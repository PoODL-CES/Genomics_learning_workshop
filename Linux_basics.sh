##### Basics
## to view file in linux
cat test_table1.txt
#or
less test_table1.txt #it will show the file contents in a new window, press q to return back to the original window

#to list contents in a directory
ls <folder_name> 

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
## 1
mkdir genomics_workshop_PoODL #creates new directory or folder named "genomics_workshop_PoODL"
mv test_table1.txt genomics_workshop_PoODL # moves the file "test_table1.txt" into the newly created directory/folder "genomics_workshop_PoODL"
cd genomics_workshop_PoODL ## changes the directory to genomics_workshop_PoODL

#### Count the number of lines in test_table1.txt

### Solution
## 1
wc -l test_table1.txt #(wc - wordcount; l- lines)

#### Count the number of characters in test_table1.txt

### Solution
## 1
wc -m test_table1.txt (wc: word count command, -m: counts characters)

#### Count the number of columns per row in test_table1.txt

### Solution
## 1
awk -F'       ' '{print NF}' test_table1.txt ### -F is the option used to declare the field separator i.e. tab here and NF is built-in variable in awk that holds the number of fields (columns) in the current line.

##2
awk '{print NF}' test_table1.txt

#3
while read line; do echo "$line" | wc -w; done < test_table1.txt

#(read: reads one line from the input and assigns it to the variable "line". the loop continues until the end.) 
#(echo "$line": This command prints the current line stored in the variable line)
#(wc -w: counts the number of columns in the current line)


#### Count the number of times "1" appears in test_table1.txt


#### Count the number of words containing "1" in test_table1.txt
