package EPSAnalysis;

public class TrackDate{

	private int year, month, date, time;

	public TrackDate(int year, int month, int date, int time){
	
		this.year = year;
		this.month = month;
		this.date = date;
		this.time = time;
	}
	
	public TrackDate(TrackDate trackdate){
		
		this.year = trackdate.getYear();
		this.month = trackdate.getMonth();
		this.date = trackdate.getDate();
		this.time = trackdate.getTime();
	}
	
	public TrackDate(String trackdate){
		
		if(trackdate.equals("0 0.000000")){
			year = 0;
			month = 0;
			date = 0;
			time = 0;

		}
		else{
			year = Integer.parseInt(trackdate.substring(0,4));
			month = Integer.parseInt(trackdate.substring(4,6));
			date = Integer.parseInt(trackdate.substring(6,8));
			time = Integer.parseInt(trackdate.substring(8,10));
		}
	}
	
	public int getYear(){	
		return year;
	}
	
	public int getMonth(){
		return month;
	}

	public int getDate(){
		return date;
	}
	
	public int getTime(){
		return time;
	}
	
	
	public void setYear(int year){
		this.year = year;
	}
	
	public void setMonth(int month){
		this.month = month;
	}
	
	public void setDate(int date){
		this.date = date;
	}
	
	public void setTime(int time){
		this.time = time;
	}
	
	public String toString(){
		String stringMonth = null;
		String stringDate = null;
		String stringTime = null;
		if(month < 10){
			stringMonth = "0" + Integer.toString(month);
		}
		else{
			stringMonth = Integer.toString(month);
		}
		if(date < 10){
			stringDate = "0" + Integer.toString(date);
		}
		else{
			stringDate = Integer.toString(date);
		}
		if(time < 10){
			stringTime = "0" + Integer.toString(time);
		}
		else{
			stringTime = Integer.toString(time);
		}
		
		return (Integer.toString(year) + stringMonth + stringDate + stringTime);
	}
	
	
	public void AddTimeFrames(int frameNo){
		
	
		for(int i=0; i< frameNo; i++){
		
		
			if((date == noOfDaysInMonth(year, month)) && (time == 18)){
			
				date = 1;
				time = 0;
				if(month == 12){
					year += 1;
					month = 1;
				}
				else{
					month +=1;
				}
			}
			else{
			
				if(time == 18){
					time = 0;
					date +=1;
			
				}
				else{
					time += 6;
				
				}
			}
		}
			
	}

	public void SubtractTimeFrames(int frameNo){
		
	
		for(int i=0; i< frameNo; i++){
			
			if((date ==1) && (time == 0)){

				time = 18;
				if(month == 1){
					date = 31;
					year -= 1;
					month = 	12;
				}
				else{
					date = noOfDaysInMonth(year, month-1);
					month -= 1;
				}
				
			}
			
			else{
			
				if(time == 0){
					time = 18;
					date -=1;
				}
				else{
					time -= 6;
				}
			}

		}
		
	}
	
	public static boolean GreaterThan(TrackDate trackdate1, TrackDate trackdate2){//returns true if trackdate1 > trackdate2, false otherwise
	
		boolean flag = true;
		if (Equal(trackdate1, trackdate2)){
			flag = false;
		}
		if(trackdate2.getYear() > trackdate1.getYear()){
			flag = false;
		}
		if(trackdate1.getYear() == trackdate2.getYear()){
				
			if(trackdate2.getMonth() > trackdate1.getMonth()){
				flag = false;
			}
			if(trackdate2.getMonth() == trackdate1.getMonth()){

				if(trackdate2.getDate() > trackdate1.getDate()){
					flag = false;
				}
				if(trackdate2.getDate() == trackdate1.getDate()){	
						
					if(trackdate2.getTime() > trackdate1.getTime()){
							flag = false;
					}
				}
			}
		}
		return flag;
		}
						
	
	public static boolean Equal(TrackDate trackdate1, TrackDate trackdate2){//returns true if trackdate1=trackdate2, false otherwise
	
		boolean flag = true;
		if(trackdate1.getYear() != trackdate2.getYear()){
			flag = false;
		}
		if(trackdate1.getMonth() != trackdate2.getMonth()){
			flag = false;
		}
		if(trackdate1.getDate() != trackdate2.getDate()){
			flag = false;
		}
		if(trackdate1.getTime() != trackdate2.getTime()){
			flag = false;
		}
		return flag;
	}

	public static int noOfDaysInMonth(int year, int month){ //returns the number of days in a given month of a given year

		switch(month) {
			
			case 1:
				return 31;
			case 2:
				if(year%4 == 0){
					return 29;
				}
				else{
					return 28;
				}
			case 3:
				return 31;
			case 4:
				return 30;
			case 5:
				return 31;
			case 6:
				return 30;
			case 7:
				return 31;
			case 8:
				return 31;
			case 9:
				return 30;
			case 10:
				return 31;
			case 11:
				return 30;
			case 12:
				return 31;
			default:
				System.err.println("Selection error " + month );
				System.exit(1);
				return -1;
		}
	}
	
	public static int leadTime(TrackDate date, TrackDate forecastStart){
		
		int leadTime = 0;
		if(GreaterThan(date, forecastStart)){
			boolean flag = true;
			while(flag){
				if(Equal(date, forecastStart)){
					flag=false;
				}
				else{
					leadTime++;
					forecastStart.AddTimeFrames(1);
				}
			}
		}
		return leadTime;

	}	

}
