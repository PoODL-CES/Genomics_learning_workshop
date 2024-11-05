##### Task 1
#### Create the following table in a file named test_table1.txt on terminal
1  2  3  4  5
6  7  8  9  10 
11  12  13  14
15  16  17  18

### Solution 
## 1
touch test_table.txt #creates a new text file in the directory
{
    printf "%2d  %2d  %2d  %2d  %2d\n" 1 2 3 4 5
    printf "%2d  %2d  %2d  %2d  %2d\n" 6 7 8 9 10
    printf "%2d  %2d  %2d  %2d\n" 11 12 13 14
    printf "%2d  %2d  %2d  %2d\n" 15 16 17 18
} > test_table.txt # ">" forwards the entered table data into test_table.txt

## 2
cat > test_table1.txt # Creates a new file in the directory or outputs the contents of the file to the terminal.
{
 1    2    3    4    5
 6    7    8    9    10
 11    12    13    14
 15    16    17    18
 }
 
## 3
nano test_table1.txt ## nano is a text editor. it will create and open a file named test_table1.txt. type the text in it. Exit with ctrl X and save with y.

#### Create a directory named genomics_workshop_PoODL and move the file test_table1.txt into it
### Solution
