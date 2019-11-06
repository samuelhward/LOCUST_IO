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
    edits prec_mod.f90 within LOCUST folder by default
---
"""


###################################################################################################
#Preamble

try:
    import sys
    import pathlib
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

if __name__=='__main__':

    try:
        import context
    except:
        raise ImportError("ERROR: context.py could not be imported!\nreturning\n")
        sys.exit(1)

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/src/support.py could not be imported!\nreturning\n") 
    sys.exit(1)

##################################################################
#Main

def LOCUST_edit_var(filepath_in='prec_mod.f90',filepath_out='prec_mod_edited.f90',**variables):
    """
    function to change hardcoded parameters initialised within a FORTRAN source file

    notes:
        only retains comments if made on first line that variable is declared (in case of spillage over multiple lines)
        formatting of edited lines will differ slightly due to blank space not being reinserted
    args:
        filepath_in - location of source file 
        filepath_out - location of output file
        variables - set of kwargs denoting variable names and values to set them to
    usage:
        python LOCUST_edit_var.py --vars c file_tet --vals 5 "'some string formatted like this'" --filepath_in prec_mod.f90
        LOCUST_edit_var(some_string="'I <3 LOCUST'") #within python - remember strings are interpreted literally so to inserting the text 'text' requires a string "'text'"
    """

    with open(filepath_in,'r') as file_input:
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
            while (True):
                if not lines[spilled_line].split(): #if line is blank then ignore
                    pass
                elif ('::' in lines[spilled_line] #delete lines until one of these characters is found - denoting a comment line, pragma line etc.
                or '#'==lines[spilled_line].split()[0][0] 
                or '!'==lines[spilled_line].split()[0][0]
                or 'stop' in lines[spilled_line]
                or 'end module' in lines[spilled_line]
                or 'end program' in lines[spilled_line]
                or 'end subroutine' in lines[spilled_line]): 
                    break
                else:
                    del(lines[spilled_line])

        with open(filepath_out,'w') as file_output: #dump results
            for line in lines:
                file_output.write(line)

if __name__=='__main__':

    try:
        import argparse
    except:
        raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
        sys.exit(1)

    parser=argparse.ArgumentParser(description='function to edit variables in LOCUST source files')
    parser.add_argument('--filepath_in',type=str,action='store',default=support.dir_locust / 'prec_mod.f90',dest='filepath_in',help="location of source file",required=False)
    parser.add_argument('--filepath_out',type=str,action='store',default=support.dir_locust / 'prec_mod_edited.f90',dest='filepath_out',help="location of output file",required=False)
    parser.add_argument('--vars',nargs='+',type=str,action='store',dest='variables',help="variables to replace in file e.g. --vars threadsPerBlock file_tet",required=True)
    parser.add_argument('--vals',nargs='+',type=str,action='store',dest='values',help='values to replace in file e.g. --vals 64 \\\'"some string value"\\\'',required=True)
    args=parser.parse_args()
    variables=dict(zip(args.variables,args.values))

    LOCUST_edit_var(filepath_in=args.filepath_in,filepath_out=args.filepath_out,**variables)

#################################

##################################################################

###################################################################################################