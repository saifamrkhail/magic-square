# The magic square challange - Group 01


> Kranzl Manuel, ai22m038 | Rahmani Saifur Rahman, ai22m055

## Run the application
This application was developped and tested on a Linux and Windows. Run tests for magic squares from 3 to 9 with `make test`.

```sh
$ make test
```

For more details on `msfinder` print help.

```sh
$ ./msfinder --help
```

The result is printed to the console and written into `magic-sqaures.txt`.
We recommend a population of 10,000, which is also the default population. A very big population slows down the programm.

## Results 
Our result are written down in the PDF. There you will find the list of the magic squares we found wiht our `msfinder`. Additionally we noted which crossover worked best for us. The strategy of elitism helped to find magic squares with certainty. Then we used  quick-sort and OpenMP to speed the algorithm. 


