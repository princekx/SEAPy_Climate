package EPSAnalysis;

public class Round{


	public static double round(double val, int places) {
	
		long factor = (long)Math.pow(10,places);

		// Shift the decimal the correct number of places
		// to the right.
		val = val * factor;

		// Round to the nearest integer.
		long tmp = Math.round(val);

		// Shift the decimal the correct number of places
		// back to the left.
		return (double)tmp / factor;
   }
	
	
	public static String round2(double val, int places){
		
		long factor = (long)Math.pow(10,places);

		val = val * factor;

		long tmp = Math.round(val);

		val = (double)tmp / factor;
		

		String number = "" + val;
		
		String test = number.substring(number.indexOf(".") +1);
		
		int a = places-test.length();
		
		if(a!=0){
		
			for(int i=0; i<a; i++){
		
				number += "0";
		
			}
		}	
		return number;
	
		}	
	


}
