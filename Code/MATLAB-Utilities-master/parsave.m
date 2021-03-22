function parsave(filename,x,varname)
eval([varname '= x;']);  
save(filename,varname)
end