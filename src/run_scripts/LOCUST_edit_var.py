#LOCUST_edit_var.py
 
"""
Samuel Ward
27/09/2019
----
edit values of initialised variables in LOCUST source code e.g. prec_mod.f90 
---
usage:
    see README.md for usage
 
notes:         
---
"""


###################################################################################################
#Preamble

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/support.py could not be imported!\nreturning\n") 
    sys.exit(1)

##################################################################
#Main

def LOCUST_edit_var(filename_in='prec_mod.f90',filename_out='prec_mod_edited.f90',**variables):
    """
    function to change initialised values within LOCUST source files

    notes:
        only retains comments if made on first line that variable is declared (in case of spillage over multiple lines)
        formatting of edited lines will differ slightly due to blank space not being reinserted
    args:
        filename_in - location of source file 
        filename_out - location of output file
        variables - set of kwargs denoting variable names and values to set them to
    usage:
        python LOCUST_edit_var.py --vars c file_tet --vals 5 \'"some string formatted like this"\' --filename_in prec_mod.f90

    """

    filename_in=support.dir_locust / filename_in #default to LOCUST directory so user can just supply filenames
    filename_out=support.dir_locust / filename_out

    with open(filename_in,'r') as file_input:
        lines=file_input.readlines()

        #interate through each line and check variable names
        spilled_lines=[] #record lines holding edited variables which spill over (&)
        for counter,line in enumerate(lines):

            if line[0] is not '#' and line[0] is not '!': #ignore preprocessor/comment lines
                if '=' in line and '::' in line:
                    line_split_equals=line.split('=')
                    number_equals_signs=len(line_split_equals)-1 #number of equals signs on this line - matters because strings sometimes have length specified with =
                    variable_name_prec_mod=line_split_equals[number_equals_signs-1].split()[-1] #holds current name of variable stored on this line in prec_mod.f90

                    for variable_name,variable_value in variables.items():
                        if variable_name == variable_name_prec_mod:

                            commented=False
                            if '!' in line_split_equals[-1]: #line is commented
                                commented=True
                                line_comment=line_split_equals[-1].split('!')[1] #extract comment on this line                    
                                line_split_equals[-1]=line_split_equals[-1].split('!')[0] #remove comment

                            if '&' in line_split_equals[-1]: #check for spill onto additional lines
                                line_split_equals[-1]=line_split_equals[1].replace('&',' ',3)
                                spilled_lines.append(counter+1)
                            
                            #create new line here
                            new_line=line_split_equals[:-1] #include all text up to before equals sign in 'var = value'
                            new_line.append(str(variable_value))
                            new_line='='.join(new_line)
                            if commented: #re-add comment if previously-existed
                                new_line=' !'.join([new_line,line_comment])
                            lines[counter]=new_line.strip()+'\n'

                            break

        for spilled_line in sorted(spilled_lines,reverse=True): #go back and delete spilled lines where we replaced variables
            while ('::' not in lines[spilled_line] #delete lines until one of these characters is found - denoting a comment line, pragma line etc.
                and '#'!=lines[spilled_line][0] 
                and '!'!=lines[spilled_line][0]
                and 'end module prec_mod' not in lines[spilled_line]): 
                del(lines[spilled_line])

        with open(filename_out,'w') as file_output: #dump results
            for line in lines:
                file_output.write(line)

if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description='function to edit prec_mod.f90 file in LOCUST source code')
    parser.add_argument('--filename_in',type=str,action='store',default='prec_mod.f90',dest='filename_in',help="location of source file",required=False)
    parser.add_argument('--filename_out',type=str,action='store',default='prec_mod_edited.f90',dest='filename_out',help="location of output file",required=False)
    parser.add_argument('--vars',nargs='+',type=str,action='store',dest='variables',help="variables to replace in prec_mod.f90 e.g. --vars threadsPerBlock file_tet",required=True)
    parser.add_argument('--vals',nargs='+',type=str,action='store',dest='values',help='values to replace in prec_mod.f90 e.g. --vals 64 \\\'"some string value"\\\'',required=True)
    args=parser.parse_args()
    variables=dict(zip(args.variables,args.values))

    LOCUST_edit_var(filename_in=args.filename_in,filename_out=args.filename_out,**variables)

#################################

##################################################################

###################################################################################################