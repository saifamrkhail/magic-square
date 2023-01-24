# The magic square challange - Group 01


> Kranzl Manuel, ai22m038 | Rahmani Saifur Rahman, ai22m055

## Run the application
This application was developped and tested on a Linux and Windows.

```sh
$ make test
```

For more details you can print help.

```sh
$ msfinder --help
```

The result is printed to the console and written into `magic-sqaures.txt`.

We recommend a population of 100,000, which is also the default population. A very big population slows down the programm.

## Results 
In the PDF magic squares are listeed again. Additionally we wrote down which crossover worked best for us.

The strategy of elitism and writing a quick sort algorithm and using OpenMP helped us speed up the algorithm. 
### Magic Sqaure 3x3
```
9 4 2 
1 8 6 
5 3 7 
```

### Magic Square 4x4
```
16 4 11 3 
9 13 10 2 
8 5 6 15 
1 12 7 14 
```

### Magic Square 5x5
```
14 20 5 4 22 
21 8 6 19 11 
12 24 13 9 7 
1 10 16 15 23 
17 3 25 18 2 
```
### Magic Square 6x6
```
33 16 21 32 5 4 
3 9 29 17 22 31 
20 8 35 27 7 14 
24 18 12 2 19 36 
25 34 1 10 30 11 
6 26 13 23 28 15 
```

### Magic Square 8x8
```
40 5 46 1 16 60 36 56 
61 18 4 58 43 8 54 14 
22 41 11 44 34 45 12 51 
7 53 49 55 32 30 15 19 
50 25 23 9 38 35 47 33 
20 28 26 29 39 17 59 42 
3 63 37 2 48 52 31 24 
57 27 64 62 10 13 6 21 
```

### Magic Square 9x9
```
40 50 51 42 74 64 26 6 16 
81 63 33 19 29 61 67 14 2 
18 41 39 31 28 9 68 56 79 
36 43 35 58 1 57 76 38 25 
7 78 66 34 30 24 27 59 44 
69 13 10 53 60 55 37 23 49 
45 3 12 80 22 75 15 52 65 
62 70 46 47 54 4 21 48 17 
11 8 77 5 71 20 32 73 72 
```

### Magic Square 10x10
```
59 9 4 90 71 96 48 53 42 33 
79 99 16 67 32 35 52 49 22 54 
15 5 87 19 81 20 72 65 43 98 
14 8 47 24 85 100 93 46 61 27 
55 97 83 62 3 10 84 73 13 25 
80 92 78 51 12 11 2 76 74 29 
64 21 36 58 91 60 7 6 94 68 
40 95 38 28 45 18 34 30 89 88 
17 56 39 75 41 86 63 70 1 57 
82 23 77 31 44 69 50 37 66 26 
```

### Magic Square 11x11
```
91 3 32 101 22 51 6 99 93 94 79 
41 90 25 23 102 29 98 42 9 105 107 
58 56 108 17 67 70 62 109 30 74 20 
106 35 12 45 77 34 117 27 43 54 121 
103 68 57 55 115 100 28 71 50 13 11 
48 104 31 76 37 52 61 86 40 72 64 
19 66 111 21 63 87 10 33 92 89 80 
8 83 49 75 110 59 116 97 26 4 44 
47 114 73 95 7 118 39 1 88 5 84 
38 36 60 85 2 18 120 82 119 65 46 
112 16 113 78 69 53 14 24 81 96 15 
```
