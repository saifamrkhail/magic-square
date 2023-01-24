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

### Magic Square 7x7
```
36 28 22 23 20 5 41
39 21 4 8 27 43 33
47 26 2 32 24 34 10
18 48 35 15 11 17 31
7 38 40 3 12 46 29
19 1 30 49 37 14 25
9 13 42 45 44 16 6
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

### Magic Square 12x12
```
29 136 61 41 79 69 51 70 82 31 96 125 
112 33 36 113 6 122 92 104 15 120 100 17 
111 133 32 80 86 44 71 18 83 68 55 89 
90 143 126 108 50 12 9 45 67 114 5 101 
97 91 76 129 66 46 14 53 141 84 60 13 
99 49 105 7 65 10 134 117 98 43 123 20 
30 87 39 42 131 140 132 52 40 19 37 121 
78 109 116 26 124 94 56 57 27 77 4 102 
130 11 75 88 28 35 2 110 119 21 144 107 
22 3 8 115 139 106 93 48 73 103 137 23 
34 59 138 74 95 64 81 142 62 72 24 25 
38 16 58 47 1 128 135 54 63 118 85 127 
```

### Magic Square 13x13
```
87 101 100 51 58 79 138 12 113 73 22 132 139 
80 55 13 64 156 65 90 36 123 135 82 122 84 
47 20 149 109 27 75 52 151 115 88 93 169 10 
69 91 114 130 28 158 99 38 86 78 67 39 108 
77 146 23 118 102 120 1 89 59 43 97 76 154 
24 127 112 148 141 3 34 32 125 145 92 16 106 
166 21 70 129 81 161 153 14 6 37 165 45 57 
30 95 9 26 136 134 164 144 98 11 107 41 110 
15 162 4 96 66 150 8 116 119 105 128 111 25 
159 17 152 7 155 19 2 157 46 103 18 133 137 
83 40 142 54 61 49 167 72 31 160 35 163 48 
147 62 143 56 44 63 126 104 60 85 68 53 94 
121 168 74 117 50 29 71 140 124 42 131 5 33 
```

### Magic Square 14x14
```
41 194 26 161 147 77 142 34 129 63 171 89 68 37 
45 108 96 84 166 114 47 25 191 141 138 21 67 136 
160 46 44 195 35 88 149 152 31 76 90 20 143 150 
169 99 176 163 133 50 72 22 55 104 3 79 122 132 
193 43 4 71 118 74 144 32 173 116 58 188 39 126 
110 145 112 115 52 78 28 155 33 156 153 106 29 107 
56 18 69 36 57 64 127 170 98 158 54 196 125 151 
140 40 128 124 1 187 102 183 81 51 9 19 134 180 
82 131 83 48 137 189 61 100 65 5 167 184 85 42 
6 172 12 70 159 121 120 139 8 111 165 16 175 105 
91 113 177 103 10 164 185 17 179 80 7 23 192 38 
117 87 181 62 186 14 13 95 135 86 190 162 2 49 
15 130 93 24 119 148 97 109 174 75 73 94 168 60 
154 53 178 123 59 11 92 146 27 157 101 182 30 66 
```

### Magic Square 15x15
```
220 8 187 39 162 95 38 171 41 91 193 215 133 18 84 
17 168 172 209 131 145 188 66 25 110 62 76 71 198 57 
123 3 67 150 181 173 92 42 59 206 2 117 203 158 119 
191 157 58 103 69 130 35 139 208 185 124 48 182 55 11 
16 142 166 13 222 81 132 4 126 121 200 51 96 118 207 
19 204 54 201 14 26 56 214 180 53 90 120 178 159 127 
30 179 213 40 65 221 224 93 137 29 15 129 34 211 75 
189 61 196 114 85 98 70 135 1 202 87 102 52 143 160 
169 28 148 197 23 10 113 174 115 21 184 24 176 88 225 
107 74 72 146 194 104 170 47 195 153 37 128 108 140 20 
100 101 12 163 217 83 199 175 94 5 64 218 80 138 46 
177 73 63 155 116 186 32 161 97 60 156 134 109 9 167 
50 86 33 77 78 79 106 216 99 205 144 183 141 149 49 
165 219 210 6 111 112 151 22 154 31 147 105 7 43 212 
122 192 44 82 27 152 89 36 164 223 190 45 125 68 136 
```