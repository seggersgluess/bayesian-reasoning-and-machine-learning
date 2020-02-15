package DataManagement;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;

public class TextDataManager {

    public String fileName;	
    public String text;
	
	public void textReader(String fileName) throws Exception {
		
		this.fileName = fileName;
		
		BufferedReader reader = getBufferedReader();
		String st;
		
		if((st = reader.readLine()) != null){
			text = Utilities.Utilities.trim(st);
			
			while((st = reader.readLine()) != null){
				if(st.isEmpty() == false) {
					text += " " + Utilities.Utilities.trim(st);
				}
			}
			
		}else {
			System.out.println("No text found in loaded file.");
		}
		
		reader.close();
	}
	
	
    public BufferedReader getBufferedReader(){
    	
		BufferedReader reader = null;
		
		try {
			reader = new BufferedReader(new FileReader(fileName));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
    	
		return reader;
		
    }
	
    
    public String getText() {
    	if(text == null) {
    		System.out.println("No text sting loaded yet.");
    		return text;
    	}else {
    		return text;
    	}
    }
    
    
    public String getfileName() {
    	return fileName;
    }
    
    
    public void printTextString() {
    	if(text == null) {
    		System.out.println("No text sting loaded yet.");
    	}else {
    		System.out.println(text);
    	}
    }
    
    
    public static void test() {
    	
    	String dirName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/NaiveBayesClassifiers/FederalistPapers/Federalist_No1.txt";
    	
    	TextDataManager textManager = new TextDataManager();
    	try {
			textManager.textReader(dirName);
		} catch (Exception e) {
			System.out.println(dirName + " not a valid file name.");
			e.printStackTrace();
		}

    	textManager.printTextString();
    }
    
	
    public static void main(String[] args) throws Exception {
    	test();    
    }
    
}
