package TextMining;

import java.util.ArrayList;

import DataManagement.TextDataManager;

public class TextDataFactory {

	String [] signs2Remove = {",", ".", ";", ":", "(", ")","[","]","{","}","|"};
	
	ArrayList<String> textAsList;
	int wordCount = 0;
	
	public String loadTextDataAsString(String dirName) {
		
		TextDataManager textManager = new TextDataManager();
		try {
			textManager.textReader(dirName);
		} catch (Exception e) {
			System.out.println(dirName + " not a valid file name.");
			e.printStackTrace();
		}
		
		String str = textManager.getText();
		return str;
	}
	
	
	public String removeSpecialSignsFromString(String str) {
		
		String cleanStr = "";
		
		int strLength = str.length()-1;	
		int blankCount = 0;
		
		for(int i=0; i<strLength; i++) {
			String s = str.substring(i, i+1);
			int [] idx = Utilities.Utilities.get_idx(signs2Remove, s);
			if(idx[0] ==-1) {	
				if(s.contentEquals(" ")==true) {
					blankCount++;
				}else {
					blankCount = 0;
				}
				if(blankCount<=1) {
					cleanStr += s;
				}
			}
		}
		
		return cleanStr;
	}
	
	
	public void getVectorizedText(String str, boolean clean) {
		
		if(clean == true) {
			str = removeSpecialSignsFromString(str);
		}
		
		textAsList = new ArrayList<String>();
		
		int strLength = str.length();
		boolean newWord = true;
		String word = null;
		
		for(int i=0; i<strLength; i++) {
			String s = str.substring(i,i+1);
			if(s.contentEquals(" ")!=true){
				if(newWord == true) {
					word = s.toLowerCase();
					newWord = false;
				}else {
					word += s.toLowerCase();
				}
			}else {
				if(newWord==false) {
					textAsList.add(word);
				}		
				newWord = true;
			}
		}
		textAsList.add(word);
		wordCount = textAsList.size();
	}
	
	
	public ArrayList<String> getTextList() {
		return textAsList;
	}
	
	
	public int getWordCount() {
		return wordCount;
	}
		
}
