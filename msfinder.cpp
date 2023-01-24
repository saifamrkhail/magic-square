#include <getopt.h>
#include <omp.h>
#include <string.h>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <random>

using namespace std;

void quickSort(int *fitness, int population, int square, int **squares, int left, int right);

class MagicSquares {
   private:
    int m_N;
    int m_square;
    int m_population;
    int **m_squares;
    int *m_fitness;

   public:
    MagicSquares(int N, int size);
    void calFit();
    void printFirst();
    int getFit(int pos);
    void breed();
    void mutate(int i);
    void copy(int i, int j);
    void crossOver(int n, int i, int j);
    void print(int i);
};

int main(int argc, char *argv[]) {
    int N = 3;
    int population = 100000;
    int del = 1000;

    /*###########################
    parse command line arguments
    ############################*/

    int arguments;
    while (true) {
        static struct option long_options[] = {
            {"dimension", required_argument, 0, 'n'},
            {"population", required_argument, 0, 'p'},
            {"help", no_argument, 0, 'h'}};
        /* getopt_long stores the option index here. */
        int option_index = 0;

        arguments =
            getopt_long(argc, argv, "hn:p:", long_options, &option_index);

        /* Detect the end of the options. */
        if (arguments == -1) {
            break;
        }

        switch (arguments) {
            case 'h': {
                cout << "usage:" << endl;
                cout << "msfinder [--dimension | -n <number>] [--population | "
                        "-p <number>]"
                     << endl;
                cout << "[--dimension | -d] <number> (the dimension of the "
                        "maqic square to search for (3x3, 4x4, 5x5, ..)), by "
                        "default 3"
                     << endl;
                cout << "[--population | -p] <number> (the population size, by "
                        "default 100000)"
                     << endl;
                exit(1);
            }

            case 'n': {
                if (strcmp("-p", optarg) == 0 ||
                    strcmp("--population", optarg) == 0) {
                    std::cerr << "Invalid argument for --dimension: " << optarg
                              << std::endl;
                    std::cerr << optarg << " is an argument of msfinder"
                              << std::endl;
                    exit(0);
                }
                try {
                    N = stoi(optarg);
                    if (N < 3 || N > 10) {
                        cout << "WARNING: N is out of range, programm might "
                                "not work."
                             << endl;
                    }
                    //adjust population size to dimension size
                    if (N == 3) {
                        population = 500;
                    } else if (N == 4) {
                        population = 1000;
                    } else if (N == 5) {
                        population = 10000;
                    }

                    if (N > 4) {
                        del = 10;
                    }
         
                } catch (const std::invalid_argument &e) {
                    std::cerr << "Invalid argument for dimension: " << optarg
                              << std::endl;
                    exit(0);
                }
                break;
            }

            case 'p': {
                if (strcmp("-n", optarg) == 0 ||
                    strcmp("--dimension", optarg) == 0) {
                    std::cerr << "Invalid argument for --population: " << optarg
                              << std::endl;
                    std::cerr << optarg << " is an argument of msfinder"
                              << std::endl;
                    exit(0);
                }
                try {
                    population = stoi(optarg);
                } catch (const std::invalid_argument &e) {
                    std::cerr << "Invalid argument for population: " << optarg
                              << std::endl;
                    exit(0);
                }
                break;
            }

            case '?': {
                for (int i = optind; i < argc + 1; i++) {
                    if (strcmp("-n", argv[i - 1]) == 0 ||
                        strcmp("--dimension", argv[i - 1]) == 0) {
                        cerr << "option --dimension requires an argument"
                             << endl;
                        exit(0);
                    } else if (strcmp("-p", argv[i - 1]) == 0 ||
                               strcmp("--population", argv[i - 1]) == 0) {
                        cerr << "option --population requires an argument"
                             << endl;
                        exit(0);
                    }
                }
                break;
            }
            default: {
                abort();
            }
        }
    }

    /* Print any remaining command line arguments (not options). */
    if (optind < argc) {
        printf("non-option argument: ");
        while (optind < argc) {
            printf("%s ", argv[optind++]);
        }
        putchar('\n');
        exit(0);
    }

    // Initialize a population of magic squares with randomly generated values
    MagicSquares squares(N, population);

    // loop until magic square is found
    long int i = 0;
    while (true) {
        // Evaluate the fitness of each square in the population,
        // using a function (sum of absolute errors of rows and columns) that
        // measures how close the square is to being a magic square
        squares.calFit();

        // Select the most fit individuals (25%) from the current population to
        // serve as parents for the next generation and breed the selected
        // parents to create a new population of magic squares, using crossover
        // and mutation to introduce genetic variation. Mutate the parents.
        squares.breed();

        // every 100 loops print the fitness of the best square for better
        // tracking of the programm
        if (i % del == 0) {
            cout << i << ":  fitness: " << squares.getFit(0) << endl;
        }

        ++i;
    }
    return 0;
}

MagicSquares::MagicSquares(int N, int size) {  // creating random squares
    m_N = N;
    m_population = size;
    m_square = m_N * m_N;
    m_squares = new int *[m_population];
    m_fitness = new int[m_population];
    random_device rd;
    mt19937 gen(rd());
    //multithreading is unsafe here 
    //#pragma omp parallel for
    for (int i = 0; i < m_population; i++) {
        m_squares[i] = new int[m_square];
        for (int j = 0; j < m_square; j++) {
            m_squares[i][j] = j + 1;
        }
        shuffle(m_squares[i], m_squares[i] + m_square, gen);
    }
}

void MagicSquares::calFit() {  // calculate the fitness of all squares
    int error = 0;
    const int magic_constant = (m_N * (m_N * m_N + 1)) / 2;
    int val;
    #pragma omp parallel for reduction(+ : error)
    for (int i = 0; i < m_population; ++i) {
        error = 0;
        // error in rows
        for (int j = 1; j < m_N; ++j) {
            val = 0;
            #pragma omp parallel for reduction(+ : val)
            for (int k = 0; k < m_N; ++k) {
                val += m_squares[i][j * m_N + k];
            }
            error += abs(magic_constant - val);
        }
        // error in columns
        for (int j = 0; j < m_N; ++j) {
            val = 0;
            #pragma omp parallel for reduction(+ : val)
            for (int k = 0; k < m_N; ++k) {
                val += m_squares[i][k * m_N + j];
            }
            error += abs(magic_constant - val);
        }
        m_fitness[i] = error;

        // if error=0 the square is a magic square, print it!
        if (error == 0) {
            this->print(i);
        }
    }
}
// print the first square, for troubleshooting
void MagicSquares::printFirst() {  
    for (int j = 0; j < m_square; ++j) {
        cout << m_squares[0][j] << " ";
    }
}

int MagicSquares::getFit(int pos) { return m_fitness[pos]; }

// note: in mutate() and copy() the fitness does not change, as it will be
// calculated in the next iteration in main().
// in swap() the fitness will be swapped also, as it is used in the sorting
// subroutine

void MagicSquares::mutate(int i) {  // mutate i by swapping two random numbers
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, m_square - 1);
    int index1 = dis(gen);
    int index2 = dis(gen);
    int buffer = m_squares[i][index1];
    m_squares[i][index1] = m_squares[i][index2];
    m_squares[i][index2] = buffer;
}

void MagicSquares::copy(int i, int j) {  // copy from i to j (j gets ereased)
    for (int k = 0; k < m_square; ++k) {
        m_squares[j][k] = m_squares[i][k];
    }
}

void MagicSquares::breed() {
    // sort the squares by their fitness from best to worst
    //this->Sort();

    quickSort(m_fitness, m_population, m_square, m_squares, 0, m_population - 1);

    // fill the later (worse) 75% of the squares with either mutation
    // of "good ones" (from the better 25%) or crossover of two "good" ones
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, 1);
    uniform_int_distribution<> dis1(0, m_population / 4 - 1);
    #pragma omp parallel for
    for (int i = m_population / 4; i < m_population; ++i) {
        if (dis(gen) == 0) {  // mutation of a "good" square (first copy then mutate)
            this->copy(dis1(gen), i);
            this->mutate(i);
        } else {  // crossover of two "good" squares
            this->crossOver(i, dis1(gen), dis1(gen));
        }
    }

    // mutate the best the worse 60% of the best 25%
    #pragma omp parallel for
    for (int i = m_population / 10; i < m_population / 4; ++i) {
        this->mutate(i);
    }
}

// save crossover of i and j in n
// more detailed description of the crossover function is in the PDF
void MagicSquares::crossOver(int n, int i, int j) {  
    int invi[m_square];
    int invj[m_square];
    int invn[m_square];
    for (int k = 0; k < m_square; ++k) {
        invi[k] = 0;
        invj[k] = 0;
    }

    // calculating inversion sequences of i and j
    for (int k1 = 1; k1 <= m_square; ++k1) {
        for (int k2 = 0; m_squares[i][k2] != k1; ++k2) {
            if (m_squares[i][k2] > k1) {
                invi[k1]++;
            }
        }
        for (int k2 = 0; m_squares[j][k2] != k1; ++k2) {
            if (m_squares[j][k2] > k1) {
                invj[k1]++;
            }
        }
    }

    // generating random pivot index
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, m_square - 1);
    int pivot = dis(gen);

    // one-point crossover of the inversion sequences
    for (int k = 0; k < pivot; ++k) {
        invn[k] = invi[k];
    }
    for (int k = pivot; k < m_square; ++k) {
        invn[k] = invj[k];
    }

    // generating n with its inversion sequence
    for (int k = 0; k < m_square; ++k) {
        m_squares[n][k] = 0;
    }
    for (int k1 = 1; k1 <= m_square; ++k1) {
        for (int k2 = 0; k2 < m_square; ++k2) {
            if (m_squares[n][k2] == 0 && invn[k1 - 1] == 0) {
                m_squares[n][k2] = k1;
                break;
            } else {
                if (m_squares[n][k2] == 0) {
                    invn[k1 - 1]--;
                } else {
                    if (invn[k1 - 1] != 0) {
                        invn[k1 - 1]--;
                    }
                }
            }
        }
    }
}

void MagicSquares::print(int i) {  // print the i-th square
    cout << "Magic Square Alarm!" << endl;
    for (int k = 0; k < m_N; ++k) {
        for (int j = 0; j < m_N; ++j) {
            cout << m_squares[i][k * m_N + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    exit(0);
}

void quickSort(int *fitness, int population, int square, int **squares, int left, int right) {
    int i = left;
    int j = right;
    int pivot = fitness[(left + right) / 2];
    while (i <= j) {
        while (fitness[i] < pivot) {
            i++;
        }
        while (fitness[j] > pivot) {
            j--;
        }
        if (i <= j) {
            swap(fitness[i], fitness[j]);
            swap(squares[i], squares[j]);
            i++;
            j--;
        }
    }
    #pragma omp parallel sections
    {
        #pragma omp section
        if (left < j) {
            quickSort(fitness, population, square, squares, left, j);
        }
        #pragma omp section
        if (i < right) {
            quickSort(fitness, population, square, squares, i, right);
        }
    }
   
}