#! /usr/bin/gawk -f

BEGIN{
#print "Change old style GET_PARAMETER(..) \
#to new Style Parameters::get(...)\n";
lineComplete=1;
cmdline="";
}
 { cmdline=cmdline$0; $0=cmdline;}
 /GET_PARAMETER/ { ##find line with GET_PARAMETER
     ## check if line is complete
     ##print " 1 "cmdline;
     gsub(/\(\)/,"XRAINERXX",cmdline );
     gsub(/\(0\)/,"XRAINERXX0",cmdline );
     gsub(/\(1\)/,"XRAINERXX1",cmdline );
     gsub(/\(2\)/,"XRAINERXX2",cmdline );
     gsub(/\(3\)/,"XRAINERXX3",cmdline );
     gsub(/\(4\)/,"XRAINERXX4",cmdline );
     split(cmdline,l,"[()]");
     if(l[3]) {lineComplete=1; }
     else {lineComplete=0;}
     if(lineComplete==1){	
	 ## find argument of GET_PARAMETER l[2]	
	 split(cmdline,l,"[()]"); 
	 split(cmdline,lw,"GET_"); ## get leading whitespaces
	 ## split arguments => a[2] is label a[3]or a[4] variable
	 split(l[2],a,",");
	 ## find variable
	 checkFormat=match(a[3],"\"");
	 if(checkFormat) {val=a[4];}
	 else{val=a[3];}
	 gsub(/^[ \t]+|[ \t]+$/,"",val)## remove whitespaces 
	 where=match(val,"&")  ##check for & 
	 if (where!=0){ 
	     gsub(/\&/,"",val);
	 }else{
	     val="(&"val")";	
	 }		 
	 cmdline= lw[1]"Parameters::get("a[2]","val");"  ;
	 gsub(/XRAINERXX/,"()",cmdline );        
	 gsub(/XRAINERXX0/,"(0)",cmdline );        
	 gsub(/XRAINERXX1/,"(1)",cmdline );        
	 gsub(/XRAINERXX2/,"(2)",cmdline );        
	 gsub(/XRAINERXX3/,"(3)",cmdline );        
	 gsub(/XRAINERXX4/,"(4)",cmdline );        
     }	  
 }
{ if(lineComplete==1) {
	gsub(/ProblemInstatVec/,"ProblemInstat",cmdline)
	gsub(/ProblemVec/,"ProblemStat",cmdline)
	gsub(/FactorLaplace_SOT/,"Simple_SOT",cmdline)
	gsub(/Laplace_SOT/,"Simple_SOT",cmdline)
	print cmdline;
	cmdline=""; } }
				

