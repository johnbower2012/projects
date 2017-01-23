#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>

std::ofstream ofile;
std::ifstream ifile;

unsigned int getFileLength(std::string fileName){
	std::ifstream file;
	file.open(fileName);
	unsigned int fileLength = 0;
	while(file.get() != EOF) fileLength++;
	file.close();
	return fileLength;
}
void fetchText(char*& text, unsigned int length, std::string fileName){
	std::ifstream file;
	file.open(fileName);
	for(int i=0;i<length;i++){
		file.get(text[i]);
	}	
	file.close();
}
void writeTextFile(char*& text, unsigned int length, std::string fileName){
	std::ofstream file;
	file.open(fileName);
	for(int i=0;i<length;i++){
		file << text[i];
	}
	file.close();
}
void printDynChar(char*& text, unsigned int length){
	for(int i=0;i<length;i++){
		std::cout << text[i];
	}
}
void encodeText(char*& iftext, unsigned int fileLength, char*& keytext, unsigned int keyLength, char*& oftext){
	for(int i=0;i<fileLength;i++){
		oftext[i] = iftext[i] + keytext[(i%keyLength)];
	}
}
void decodeText(char*& iftext, unsigned int fileLength, char*& keytext, unsigned int keyLength, char*& oftext){
	for(int i=0;i<fileLength;i++){
		oftext[i] = iftext[i] - keytext[(i%keyLength)];
	}
}


int main(int argc, char *argv[]){
	
	std::string ifname, keyname, ofname;
	unsigned int code;
/*
	if(argc<4){
		std::cout << "Bad usage. Enter also 'inFileName keyFileName outFileName.'" << std::endl;
		exit(1);
	}
	else{
		ifname = argv[1];
		keyname = argv[2];
		ofname = argv[3];
	}
*/

	std::cout << "Decode (1) or Encode (0):		";
	std::cin >> code;
	if(code!=0&&code!=1){
	 	std::cout << "Improper entry. Enter 1 to decode or 0 to encode file. Terminating program..." << std::endl;
		exit(1);
	}
	std::cout << "Enter input text file:	";
	std::cin >> ifname;
	std::cout << "Enter key text file:	";
	std::cin >> keyname;
	std::cout << "Enter output text file:	";
	std::cin >> ofname;
	
	unsigned int fileLength = getFileLength(ifname);
	unsigned int keyLength = getFileLength(keyname);

	char* iftext = new char[fileLength];
	char* oftext = new char[fileLength];
	char* keytext = new char[keyLength];

	fetchText(iftext,fileLength,ifname);
	fetchText(keytext,keyLength,keyname);

	if(code==1){
		decodeText(iftext,fileLength,keytext,keyLength,oftext);
	}
	else if(code==0){
		encodeText(iftext,fileLength,keytext,keyLength,oftext);
	}

	writeTextFile(oftext,fileLength,ofname);

	delete[] iftext, oftext, keytext;
	
	return 0;
}
