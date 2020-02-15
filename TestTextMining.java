package TextMining;

import java.util.ArrayList;
import java.util.HashMap;

public class TestTextMining {

	//Test 1: Loading and vectorizing text from txt-file in TextFactory
	public static void test1() {
    	
    	String dirName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/NaiveBayesClassifiers/FederalistPapers/Federalist_No1.txt";
    	
    	TextDataFactory textFact = new TextDataFactory();
    	String str = null;
    	try {
			str = textFact.loadTextDataAsString(dirName);
		} catch (Exception e) {
			System.out.println(dirName + " not a valid file name.");
			e.printStackTrace();
		}

    	textFact.getVectorizedText(str,true);
    	
    	ArrayList<String> textList = textFact.getTextList();
    	System.out.println(textList);
    	
    }
    
	
	//Test 2: Preparing text loaded from file and make simple statistics on these text sample
	public static void test2() {
		
    	String dirName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/NaiveBayesClassifiers/FederalistPapers/Federalist_No1.txt";
    	
		TextAnalyzer textAnalyzer = new TextAnalyzer();
		
		textAnalyzer.loadAndPrepareText(dirName);
		
		String [] words = {"people", "state", "the", "is"};
		
		//Count specified words in supplied text
		HashMap<String, Integer> wc = textAnalyzer.getWordCount(words);		
		System.out.println(wc);
		
		//Count all (unique) words in supplied text
		wc = textAnalyzer.getTextStatistics();		
		System.out.println(wc);		
				
	}
	
	
    public static void main(String[] args) throws Exception {
    	test2();    
    }
	
}
