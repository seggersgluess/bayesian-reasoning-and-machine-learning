package DataManagement;

import java.util.ArrayList;
import java.util.List;
import Utilities.Utilities;

public class JsonObject {

	static String json;
	static List<String> jsonKeys              = new ArrayList<String>();
	static ArrayList<List<String>> jsonValues = new ArrayList<List<String>>();	
	
	static List<Integer> jsonValueIdxs  = new ArrayList<Integer>();
	static List<Integer> jsonArrayIdxs  = new ArrayList<Integer>();
	
	//static ArrayList<JsonObject> jsonObjects = new ArrayList<JsonObject>();
		
	@SuppressWarnings("static-access")
	public JsonObject(String json){
		this.json = json;
	}
	
	
	public static void parseJson(){
		
		int jsonValueCount = 0;
		int quoteCount = 0;
		boolean keyAlreadySet = false;
		String str     = "";
		
		for(int i=1; i<json.length()+1; i++){
			
			String token = json.substring(i-1,i);
			
			if(isQuote(token) == true){
				quoteCount++;				
			}else{
				
				if(quoteCount == 1){					
					str=str+token;						
				}
				
				if(quoteCount == 2){										
					
					if(isColon(token) == true){
						//json key
						jsonKeys.add(str);
						keyAlreadySet = true;
						str = "";
					}
					
					if(isComma(token) == true || isRightBrace(token) == true){
						//json value is string
						List<String> values = new ArrayList<String>(1);						
						values.add(str);
						
						jsonValues.add(values);						
						jsonValueIdxs.add(jsonValueCount);
						jsonValueCount++;
						keyAlreadySet = false;
						str = "";
					}
													
					quoteCount = 0;
					
				}
				
				if(isLeftBracket(token) == true){
					//json array
					i=jsonSimpleArray(i);
					jsonArrayIdxs.add(jsonValueCount);
					jsonValueCount++;
					token = json.substring(i-1,i);
					keyAlreadySet = false;
				}
				
			}		
			
			if(keyAlreadySet == true && isLeftBrace(token) == false && quoteCount == 0){
				//json value is number
				if(isComma(token) == false && isColon(token) == false && isRightBrace(token) == false){
					str=str+token;
				}
				
				if(isComma(token) == true || isRightBrace(token) == true){
					//json value
					List<String> values = new ArrayList<String>(1);					
					values.add(str);
					
					jsonValues.add(values);
					jsonValueIdxs.add(jsonValueCount);
					jsonValueCount++;
					str = "";
					keyAlreadySet = false;
				}
				
			}
			
		}
		
	}

	
	public static int jsonSimpleArray(int posIdx){
		
		String str = "";
		boolean endLoop         = false;
		boolean valueAlreadySet = false;
		int quoteCount = 0;
		
		List<String> array = new ArrayList<String>();
		
		String token = json.substring(posIdx-1,posIdx);
		
		while(endLoop == false){
						
			if(isQuote(token) == true){
				quoteCount++;				
			}else{
				
				if(quoteCount == 1){					
					str=str+token;						
				}
				
				if(quoteCount == 2){										
									
					if(isComma(token) == true || isRightBracket(token) == true){
						//json value is string
						array.add(str);	
						valueAlreadySet = true;
						str = "";
							
						if(isRightBracket(token) == true){
							endLoop = true;
						}
						
					}
													
					quoteCount = 0;
					
				}
				
			}
		
			if(valueAlreadySet == false && isLeftBracket(token) == false && quoteCount == 0){
				//json value is number
				if(isComma(token) == false && isColon(token) == false && isRightBracket(token) == false){
					str=str+token;
				}
				
				if(isComma(token) == true || isRightBracket(token) == true){
					//json value				
					array.add(str);		
					str = "";
					
					if(isRightBracket(token) == true){
						endLoop = true;
					}
					
				}
				
			}
			
			posIdx++;
			token = json.substring(posIdx-1,posIdx);
						
		}
		
		jsonValues.add(array);
		
		return posIdx;
		
	}
	
	
	public static boolean isQuote(String token){		
		return token.contentEquals("\"");		
	}
	
	
	public static boolean isColon(String token){
		return token.contentEquals(":");
	}
	
	
	public static boolean isComma(String token){
		return token.contentEquals(",");
	}
	
	
	public static boolean isLeftBrace(String token){
		return token.contentEquals("{");
	}
	
	
	public static boolean isRightBrace(String token){
		return token.contentEquals("}");
	}
	
	
	public static boolean isLeftBracket(String token){
		return token.contentEquals("[");
	}
	
	
	public static boolean isRightBracket(String token){
		return token.contentEquals("]");
	}
	
	
	public static List<String> getValues(String key){
		
		int [] idx = Utilities.get_idx(jsonKeys, key);
		
		if(idx[0] == -1){
			throw new RuntimeException(key + " no valid key of supplied JSON file.");
		}
		
		return jsonValues.get(idx[0]);
		
	}
	
	
	public static int getNumberOfArrays(){
		return jsonArrayIdxs.size();
	}
	
	
	public static int getNumberOfValues(){
		return jsonValueIdxs.size();
	}
	
	
	public static void main(String[] args){
		
		String a1 = "{\"Key1\":123,\"Key2\":\"AAA\",\"Key3\":\"BBBBBBBB\", \"Key4\":47957}"; 
		String a2 = "{\"Key1\":\"BBBBB\",\"Key2\":\"AAA\"}";
		String a3 = "{\"Key1\":\"908.01\",\"Key2\":123,\"Key3\":\"AAA\",\"Key4\":456,\"Key5\":\"AAA\",\"Key6\":\"BALUEORL:DLUR\",\"Key7\":\"31.01.2018\"}";
		String a4 = "{\"Key1\":\"AAA\"}";
		String a5 = "{\"Key1\":[\"AAA\",\"bboal\"],\"Key2\":123,}";
		
		String a6 = "{"
				    + "\"Key1\":[123,45,66],"
				    + "\"Key2\":123,"
				    + "\"Key3\":[\"A\",\"B\",\"C\"],"
				    + "\"Key4\":[\"Audi\",\"BMW\",\"Mercedes\"]"
				    + "}";
		
		json = a6;
			
		parseJson();
		
		System.out.println(getValues("Key4"));

		
	}
	
	
}
