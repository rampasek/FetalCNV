#!/usr/bin/awk -f
{
if (match($8, /^DP=/)){	
	#split($8, a, ";");
	#if (match(a[1], /^DP=/)){
	if ($6 >= qlimit){
		print $0;
	}
}
}
