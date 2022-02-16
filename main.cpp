#include <fstream>
/*
 * Created by Santiago Fernandez 9-1-21
 * Project 1:
 *
 * given a database and data file of DNA and STRs store, process and search this
 * information.
 *
 */

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "ourvector.h"
using namespace std;

// this structure will store each DNA profile with it's respectful data
struct DNAprofile {
  string name;
  ourvector<int> SequenceFrequency;
  ourvector<ourvector<char>> dnaStr;
};

// fuction will open file and check for errors
int openF(ifstream& Input) {
  string DBname;
  cin >> DBname;

  Input.open(DBname);

  if (!Input.is_open() | Input.fail()) {
    cout << "Error: unable to open '" << DBname << "'\n";
    return 1;
  } else
    return 0;
}

// sinqData is a helper function for load_db it retrives the database DNA STRs
// and stores them in ourvecto<char> variable, checks for end of file. This will
// help create my template/database of STRs
void sinqData(ifstream& Input, ourvector<ourvector<char>>& DNAvec) {
  string tempString;

  getline(Input, tempString);

  stringstream tempSS(tempString);

  // discard 'name'
  getline(tempSS, tempString, ',');

  while (tempSS.eof() != true) {
    getline(tempSS, tempString, ',');

    ourvector<char> currentDNA;
    for (size_t i = 0; i < tempString.size(); i++) {
      currentDNA.push_back(tempString[i]);
    }

    DNAvec.push_back(currentDNA);
  }

  tempSS.clear();
}

// (MS1) load-db will read a file and create DNA profiles with the data on the
// file
int load_db(ourvector<DNAprofile>& database,
            ourvector<ourvector<char>>& dnaSamples) {
  cout << "Loading database...\n";
  database.clear();
  dnaSamples.clear();

  ifstream dbFile;
  if (openF(dbFile) == 1) return 1;

  // this ourvector will store the DNA strands that are present in the file in
  // the order they are given
  ourvector<ourvector<char>> DNAvec;
  sinqData(dbFile, DNAvec);

  // I used stringstream to parse information one line at a time, each line will
  // be stored in tempString.
  string tempString;

  while (dbFile.eof() != true) {
    DNAprofile ThisProfile;
    ThisProfile.dnaStr = DNAvec;
    dnaSamples = DNAvec;

    getline(dbFile, tempString);
    stringstream tempSS(tempString);

    // from the string stream, information is then parsed by the ',' character.
    getline(tempSS, tempString, ',');
    ThisProfile.name = tempString;

    while (tempSS.eof() != true) {
      getline(tempSS, tempString, ',');

      ThisProfile.SequenceFrequency.push_back(stoi(tempString));
    }

    tempSS.clear();

    database.push_back(ThisProfile);
  }

  dbFile.close();
  return 0;
}

// (MS3) load_dna will take in data from a file and store it in an ourvector. it
// will take all characters one by one. This will be our long, un-processed DNA
// data.
int load_dna(ourvector<char>& DNA) {
  ifstream DNAfile;
  char ch;

  cout << "Loading DNA...\n";
  DNA.clear();

  // ifstream dbFile;
  if (openF(DNAfile) == 1) return 1;

  while (DNAfile.eof() != true) {
    DNAfile.get(ch);
    DNA.push_back(ch);
  }

  return 0;
}

// (MS2) overloaded display function. This function will display the database of
// profiles will iterate through the ourvector of DNA profiles and display the
// structures data
void display(ourvector<DNAprofile>& database) {
  if (!database.size()) {
    cout << "No database loaded.\n";
  } else {
    cout << "Database loaded: \n";
    for (int i = 0; i < database.size(); i++) {
      cout << database[i].name;
      for (int j = 0; j < database[i].SequenceFrequency.size(); j++) {
        cout << " " << database[i].SequenceFrequency[j];
      }
      cout << endl;
    }
  }
}

// (MS2) overloaded display function. This function iterates through the
// un-processed ourvector of DNA and displays it
void display(ourvector<char>& DNA) {
  if (!DNA.size()) {
    cout << "No DNA loaded.\n";
  } else {
    cout << "DNA loaded: \n";
    for (int y = 0; y < (DNA.size() - 1); y++) cout << DNA[y];
  }
  cout << endl;
}

// (MS2) overloaded display function. This function will iterate through the
// processed DNA ourvector and display how many consecutive STRs were in the DNA
// outvector/ file
void display(ourvector<int>& ProcessedSTR, ourvector<ourvector<char>>& STRdb) {
  if (!ProcessedSTR.size()) {
    cout << "No DNA has been processed.";
  } else {
    cout << "DNA processed, STR counts: \n";
    for (int i = 0; i < ProcessedSTR.size(); i++) {
      for (int j = 0; j < STRdb[i].size(); j++) cout << STRdb[i][j];
      cout << ": " << ProcessedSTR[i] << endl;
    }
  }
  cout << endl;
}

// helper function for the process function. Identifies if a match id found at
// the position given.
bool match(ourvector<char>& STR, ourvector<char>& DNA, int i) {
  for (int j = 0; j < STR.size(); j++) {
    if (i >= (DNA.size() - 1)) {
      return false;
    }

    if (STR[j] != DNA[i]) {
      return false;
    } else {
      i++;
    }
  }

  return true;
}

// (MS4) process function. will process take the long un-processed DNA vector
// and the 1 possible STR and return value of the longest consecutive sequence
// of STRs are found in the DNA ourvector
int process(ourvector<char>& STR, ourvector<char>& DNA,
            ourvector<ourvector<int>>& allPositions) {
  int consCounter = 0;
  int largestCons = 0;
  bool consecutives = false;
  ourvector<int> positions;

  // finds the first matching character in STR and DNA then calls the match()
  // function to confirm there is a matching STR at that position. pushes back
  // the position into a vector of integers(positions)
  for (int i = 0; i < (DNA.size() - 2); i++) {
    if (DNA[i] == STR[0]) {
      if (match(STR, DNA, i) == true) {
        largestCons = 1;
        positions.push_back(i);
        i = i + (STR.size() - 1);
      }
    }
  }
  // pushes current STRs match posintions into ouvvector of ourvector<int>
  allPositions.push_back(positions);

  // checks the vector of positions to find if they are consecutive.
  for (int j = 1; j < positions.size(); j++) {
    if (positions[j - 1] == (positions[j] - STR.size())) {
      consCounter++;
      consecutives = true;
      if (consCounter > largestCons) {
        largestCons = consCounter;
      }
    } else
      consCounter = 0;
  }

  if (consecutives == true) {
    return (largestCons + 1);
  }
  return largestCons;
}

// overloaded helper function for Process, checks if database or DNA have been
// loaded.
int checkSize(ourvector<ourvector<char>>& STR, ourvector<char>& DNA) {
  if (!STR.size()) {
    cout << "No database loaded.\n";
    return 1;
  } else if (!DNA.size()) {
    cout << "No DNA loaded.\n";
    return 1;
  } else {
    cout << "Processing DNA...\n";
  }

  return 0;
}

// overloaded helper function for Search, checks if database or DNA have been
// loaded and if DNA has been processed.
int checkSize(ourvector<DNAprofile>& database, ourvector<int>& processedSTR,
              ourvector<char>& DNA) {
  if (!database.size()) {
    cout << "No database loaded.\n";
    return 1;
  } else if (!DNA.size()) {
    cout << "No DNA loaded.\n";
    return 1;
  } else if (!processedSTR.size()) {
    cout << "No DNA processed.\n";
    return 1;
  } else {
    cout << "Searching database...\n";
  }

  return 0;
}

// helper function for Process, processes each STR to find the longest sequence
// of each STR, then loades them into an ourvector
ourvector<int> runProcess(ourvector<ourvector<char>>& STR, ourvector<char>& DNA,
                          ourvector<ourvector<int>>& allPositions) {
  ourvector<int> longestSTR;

  if (checkSize(STR, DNA) == 0) {
    for (int n = 0; n < STR.size(); n++)
      longestSTR.push_back(process(STR[n], DNA, allPositions));
  }

  return longestSTR;
}

// (MS5) helper function for search. compares each profile to the processed DNA
// ourvector and returns true if it is a match.
bool searchProfile(ourvector<int>& dbSTRFrequency, ourvector<int>& processedSTR) {
  bool match = false;

  for (int i = 0; i < dbSTRFrequency.size(); i++) {
    if (dbSTRFrequency[i] == processedSTR[i])
      match = true;
    else
      return false;
  }
  return match;
}

// calls searchProfile() and displays if match is found.
void search(ourvector<DNAprofile>& database, ourvector<int>& processedSTR,
            ourvector<char>& DNA) {
  bool matchFound;

  if (checkSize(database, processedSTR, DNA) == 0) {
    for (int i = 0; i < database.size(); i++) {
      matchFound = searchProfile(database[i].SequenceFrequency, processedSTR);

      if (matchFound == true) {
        cout << "Found in database! ";
        cout << "DNA matches: " << database[i].name << endl;
        i = database.size();
      }
    }
    if (matchFound == false) cout << "Not found in database.\n";
  }
}

//(MS6) will return all instances of each str
void filter(ourvector<ourvector<char>>& STR,
            ourvector<ourvector<int>>& allPositions, ourvector<char>& DNA) {
  if (checkSize(STR, DNA) == 0) {
    for (int i = 0; i < STR.size(); i++) {
      cout << "Instances of ";
      for (int j = 0; j < STR[i].size(); j++) {
        cout << STR[i][j];
      }
      cout << ": " << allPositions[i].size() << endl;
    }
  }
}

// menu function manages user input and calls functions depending on input.
void menu(string& command, ourvector<DNAprofile>& database,
          ourvector<ourvector<char>>& STR, ourvector<char>& DNA,
          ourvector<int>& processedSTR,
          ourvector<ourvector<int>>& allPositions) {
  while (command != "#") {
    if (command == "load_db") {
      load_db(database, STR);
      cout << "Enter command or # to exit: ";
      cin >> command;
    } else if (command == "display") {
      display(database);
      display(DNA);
      display(processedSTR, STR);
      cout << "Enter command or # to exit: ";
      cin >> command;
    } else if (command == "load_dna") {
      load_dna(DNA);
      cout << "Enter command or # to exit: ";
      cin >> command;
    } else if (command == "process") {
      processedSTR = runProcess(STR, DNA, allPositions);
      cout << "Enter command or # to exit: ";
      cin >> command;
    } else if (command == "search") {
      search(database, processedSTR, DNA);
      cout << "Enter command or # to exit: ";
      cin >> command;
    } else if (command == "filter") {
      filter(STR, allPositions, DNA);
      cout << "Enter command or # to exit: ";
      cin >> command;

    } else {
      cout << endl << "Invalid input... \n";
      cout << endl << "Enter command or # to exit: ";
      cin >> command;
      cout << endl;
    }
  }
}

// main function calls menu function and initializes some necessary ourvectors
// creative component(MS6) filter function, (called in menu) returns all
// instances of
int main() {
  string command;
  ourvector<DNAprofile> database;
  ourvector<ourvector<char>> STR;
  ourvector<ourvector<int>> allPositions;
  ourvector<char> DNA;
  ourvector<int> processedSTR;

  cout << "Welcome to the DNA Profiling Application." << endl;
  cout << "Enter command or # to exit: ";
  cin >> command;

  menu(command, database, STR, DNA, processedSTR, allPositions);

  return 0;
}