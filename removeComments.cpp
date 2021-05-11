// C++ program to remove comments from a C/C++ program 
#include <iostream> 
#include <fstream>
#include <sstream>
#include <vector>
using namespace std; 

string removeComments(string prgm) 
{ 
	int n = prgm.length(); 
	string res; 

	// Flags to indicate that single line and multpile line comments 
	// have started or not. 
	bool s_cmt = false; 
	bool m_cmt = false; 


	// Traverse the given program 
	for (int i=0; i<n; i++) 
	{ 
		// If single line comment flag is on, then check for end of it 
		if (s_cmt == true && prgm[i] == '\n') 
			s_cmt = false; 

		// If multiple line comment is on, then check for end of it 
		else if (m_cmt == true && prgm[i] == '*' && prgm[i+1] == '/') 
			m_cmt = false, i++; 

		// If this character is in a comment, ignore it 
		else if (s_cmt || m_cmt) 
			continue; 

		// Check for beginning of comments and set the approproate flags 
		else if (prgm[i] == '/' && prgm[i+1] == '/') 
			s_cmt = true, i++; 
		else if (prgm[i] == '/' && prgm[i+1] == '*') 
			m_cmt = true, i++; 

		// If current character is a non-comment character, append it to res 
		else res += prgm[i]; 
	} 
	return res; 
} 


