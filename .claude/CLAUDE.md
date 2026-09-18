# Python formatting
* Format the Python function calls and definitions, list and dictionary definitions etc. with all arguments aligned  
horizontally on a single line, rather than stacking them vertically one per line. 
* Follow PEP 8 guidelines for long argument lists: hang subsequent lines at the opening parenthesis position  
without vertical alignment, but prefer a single horizontal line if it fits within the 120-character limit.
* Define log, error, and exception message strings in its own line before passing them as arguments.
* If a dictionary needs to be passed as an argument, define it outside the function call.
* Do not leave a single opening or closing bracket in a line
* Use f-strings for string formatting and interpolation. Do not use modulo operator (%) or `.format()`. 
  When interpolating floats format to two digital places.
* Use `typing` module for type hints wherever possible


# Plan Mode Rules
When explicitly asked to create plan, or when in a plan mode, but not otherwise, read the file `.claude/plan_mode_rules.
md` 
