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
