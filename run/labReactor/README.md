## Running validations  
Prepare cases by executing corresponding `Allrun` scripts:  
`Allrun_Leion`  
`Allrun_Berguerand`  
`Allrun_density`  
They will print commands to run corresponding cases

## Running single case  
`cp -r labReactor labReactor_act_CH4_CASENAME`  
`runCase labReactor_act_CH4 --chemistry act --gas CH4`  

You can specify any entry from caseParamsDict as -- arguement  