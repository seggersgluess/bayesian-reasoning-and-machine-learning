package TextMining;

import java.util.ArrayList;
import java.util.HashMap;

public class TextAnalyzer {

	public TextDataFactory textFactory;
	
	
	public void loadAndPrepareText(String dirName) {
		
    	textFactory = new TextDataFactory();
    	String str = null;
    	try {
			str = textFactory.loadTextDataAsString(dirName);
		} catch (Exception e) {
			System.out.println(dirName + " not a valid file name.");
			e.printStackTrace();
		}

    	textFactory.getVectorizedText(str,true);
		
	}
	
	
	public HashMap<String, Integer> getWordCount(String [] words) {
		
		if(textFactory == null) {
			throw new RuntimeException("No loaded & prepared text found for counting words.");	
		}
		
		ArrayList<String> textList = textFactory.getTextList();
		
		if(textList == null) {
			throw new RuntimeException("No loaded & prepared text found for counting words.");
		}
		
		int nWords = words.length;
		int nWordsInText = textFactory.getWordCount();
		
		HashMap<String, Integer> wordCounts = new HashMap<String, Integer>(nWords);
		
		for(int i=0; i<nWords; i++) {
			words[i].toLowerCase();
			wordCounts.put(words[i], 0);
		}
		
		for(int w=0; w<nWordsInText; w++) {
			for(int i=0; i<nWords; i++) {
				if(textList.get(w).contentEquals(words[i])) {
					wordCounts.put(words[i],wordCounts.get(words[i])+1);
				}
			}
		}
		
		return wordCounts;
	}
	
	
	public HashMap<String, Integer> getWordCount(ArrayList<String> words) {
		
		if(textFactory == null) {
			throw new RuntimeException("No loaded & prepared text found for counting words.");	
		}
		
		ArrayList<String> textList = textFactory.getTextList();
		
		if(textList == null) {
			throw new RuntimeException("No loaded & prepared text found for counting words.");
		}
		
		int nWords = words.size();
		int nWordsInText = textFactory.getWordCount();
		
		HashMap<String, Integer> wordCounts = new HashMap<String, Integer>(nWords);
		
		for(int i=0; i<nWords; i++) {
			words.get(i).toLowerCase();
			wordCounts.put(words.get(i), 0);
		}
		
		for(int w=0; w<nWordsInText; w++) {
			for(int i=0; i<nWords; i++) {
				if(textList.get(w).contentEquals(words.get(i))) {
					wordCounts.put(words.get(i),wordCounts.get(words.get(i))+1);
				}
			}
		}
		
		return wordCounts;
	}
	
	
	public HashMap<String, Integer> getTextStatistics() {
		
		if(textFactory == null) {
			throw new RuntimeException("No loaded & prepared text found for counting words.");	
		}
		
		ArrayList<String> textList = textFactory.getTextList();
		
		if(textList == null) {
			throw new RuntimeException("No loaded & prepared text found for counting words.");
		}
		
		ArrayList<String> uniqueWords = Utilities.Utilities.get_unique_elements(textList);
		
		HashMap<String, Integer> wordStats = getWordCount(uniqueWords);
			
		return wordStats;
	}
	
}
