#LOCUST_edit_var.py
 
'''
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
'''


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
    import run_scripts.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/utils.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/src/support.py could not be imported!\nreturning\n") 
    sys.exit(1)

##################################################################
#Main

def LOCUST_edit_var(filepath_in='prec_mod.f90',filepath_out='prec_mod_edited.f90',**settings):
    """
    function to change hardcoded parameters initialised within a FORTRAN source file

    notes:
        only retains comments if made on first line that variable is declared (in case of spillage over multiple lines)
        formatting of edited lines will differ slightly due to blank space not being reinserted
        assumes lines may only contain a single !, signifying a comment
    args:
        filepath_in - location of source file 
        filepath_out - location of output file
        settings - set of kwargs denoting variable names and values to set them to
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

                    commented=False
                    line_comment=''
                    line_split_comment=line.strip().split('!') #extract comment on this line  
                    if not len(line_split_comment)==1: #if there is a comment on the line
                        commented=True
                        line_comment=line_split_comment[-1]
                        number_equals_signs_comments=len(line_comment.split('='))-1 #number of equals signs in comment of this line
                        number_equals_signs-=number_equals_signs_comments

                    variable_name_prec_mod=line_split_equals[number_equals_signs-1].split()[-1] #holds current name of variable stored on this line in prec_mod.f90

                    for variable_name,variable_value in settings.items():
                        if variable_name == variable_name_prec_mod:

                            #need to take everything before comment section, split according to equals and leave off the very last entry (must preserve = that come before the = we care about)                 
                            if not commented: 
                                left_side_of_equals='='.join(substring for substring in line_split_equals[:-1])    
                                right_side_of_equals=line_split_equals[-1]    
                            else:
                                left_side_of_equals='='.join(substring for substring in line_split_comment[0].split('=')[:-1])
                                right_side_of_equals=line_split_comment[0].split('=')[-1]

                            if '&' in right_side_of_equals: #check for spill onto additional lines
                                spilled_lines.append(counter+1)
                            
                            #create new line here
                            new_line=''.join([
                            left_side_of_equals,'=',str(variable_value),' ! ',line_comment,' *edited with LOCUST_IO*']) #include all text up to before equals sign in 'var = value'
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
    parser.add_argument('--settings',nargs='+',type=str,action='store',dest='settings',help="variables and corresponding values to replace in file e.g. --vars threadsPerBlock=32 file_tet=\"\'some_filename\'\"",required=True)
    
    args=parser.parse_args()
    args.settings=run_scripts.utils.command_line_arg_parse_dict(args.settings) #transform into dictionary

    LOCUST_edit_var(filepath_in=args.filepath_in,filepath_out=args.filepath_out,**args.settings)

#################################

##################################################################

###################################################################################################