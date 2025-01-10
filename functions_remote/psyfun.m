function p = psyfun(refSRT,refslope,Lvec,maxSI) 

p = maxSI*1./(1+exp(4*refslope*(refSRT-Lvec))); 

end