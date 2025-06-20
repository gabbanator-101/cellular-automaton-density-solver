#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/* #DEFINE'S -----------------------------------------------------------------*/
#define SDELIM "==STAGE %d============================\n"   // stage delimiter
#define MDELIM "-------------------------------------\n"    // delimiter of -'s
#define THEEND "==THE END============================\n"    // end message
#define SPACE_DELIMITER " " //delimiter of one space

#define ON_OFF_STATE_COUNTS "#ON=%d #OFF=%d CELL#%d START@%d\n"
#define RULE_STEPS_FORMAT "RULE: %d; STEPS: %d.\n"
#define SOLUTION_DCP_FORMAT "AT T=%d: #ON/#CELLS %c 1/2\n"

#define CRTRNC '\r'     // carriage return character
#define NEWLINE_CHAR '\n' // newline character
#define NBRHDS 8        // number of possible neighborhoods

#define NBRHD_ELEMENTS 3  // the number of elements in each neighbourhood

#define HEADER_ROW 1  // Rows for the initial starting elementary CA

#define FIRST_INDEX 0
#define SECOND_INDEX 1
#define THIRD_INDEX 2
#define LAST_INDEX -1 // the index required to be added to obtain the last index
#define FIRST_STATE_INDEX 1 // states will be stored in array from index 1 

/* The elementary CA will be stored from index 1 to size (incl.) so NEXT_INDEX
must be added to ensure correct indexing */
#define NEXT_INDEX 1 
#define PREVIOUS_INDEX -1 

#define WRAPPING 2 // Wrap array ends
#define BINARY_BASE 2

#define ONSTATE '*' // The character representing a cell which is "ON"
#define OFFSTATE '.' // The character representing "OFF"

#define ON 1 // Integer representation of "ON"
#define OFF 0 // Integer representation of "OFF"

#define SAME 1 // Used to indicate that two bits are the same for comparison
#define DIFFERENT 0 // Used to indicates that two elements are different

#define INITIAL 0 // Initial value for variables

/* operators used in the display of the solution to DCP */
#define LESS '<'
#define GREATER '>'
#define EQUIVALENT '='

#define SOLVE_DCP_RULE1 184 // rule used as step 1 to solve the DCP
#define SOLVE_DCP_RULE2 232 // rule used as step 2 to solve the DCP

/* TYPE DEFINITIONS ----------------------------------------------------------*/
typedef int cells_t;            // base type to store states of cells

// base type to store binary representation of each neighbourhood
typedef cells_t neighbourhood_t[NBRHD_ELEMENTS];  

/* base type to store each neighbourhood, with the binary representation as well
as the updated state at time t+1 for a given neighbourhood */
typedef struct {
    neighbourhood_t binary;
    cells_t updated_state;
} neighbourhood_update_t;

// base type to store the number of evolved steps of the elementary CA
typedef struct {
    int time_step; // number of times the initial elementary CA is evolved
    int n_steps_184; // number of evolved times (rule 184)
    int m_steps_232; // number of evolved times (rule 232)
} steps_t; 

/* USEFUL FUNCTIONS ----------------------------------------------------------*/

cells_t * print_stage0_output(neighbourhood_update_t *, int *);
void determine_update_states(neighbourhood_update_t *, int);
void print_rule(neighbourhood_update_t *); 
cells_t * malloc_elementaryCA (int);
void read_elementaryCA(cells_t *, int);
int my_getchar(void);
void print_elementaryCA(cells_t *, int, int);
char convert_int_to_char(int);
void wrap_array_ends(cells_t *, int);
void convert_decimal_to_binary(neighbourhood_update_t *, int);

cells_t ** print_stage1_output(neighbourhood_update_t *, cells_t *, int, int *);
void add_elementaryCA_to_array(cells_t **, int, neighbourhood_update_t *, 
                               int, int);
cells_t * compute_next_time_CA(neighbourhood_update_t *, cells_t *, int);
int match_neighbourhoods(neighbourhood_update_t *, cells_t *);
int compare_binary_bits(cells_t *, cells_t *);
void print_on_off_state_counts(cells_t **, int);

cells_t ** print_stage2_output(cells_t **, int, steps_t *);
void solve_density_problem(cells_t **, int, int, int, int);
char find_operator(cells_t *);

void free_elementaryCA(cells_t **, steps_t *);

/* WHERE IT ALL HAPPENS ------------------------------------------------------*/
int main(int argc, char *argv[]) {
    
    // define an array of 8 neighbourhoods, as defined above
    neighbourhood_update_t neighbourhoods[NBRHDS];

    int size, time_step = INITIAL;

    // initialise elementary CA with input, while keeping track of size in main
    cells_t * elementaryCA = print_stage0_output(neighbourhoods, &size);
        
    // store the 2D dynamic array of evolved elementary CAs in stage 1
    cells_t ** all_elementaryCA = print_stage1_output(neighbourhoods, 
                                        elementaryCA, size, &time_step);

    // a structure to store the number of times the CA has evolved by rule
    steps_t evolved_steps; 
    evolved_steps.time_step = time_step;

    // store the updated 2D array of elementary CAs after rules 184 and 232
    all_elementaryCA = print_stage2_output(all_elementaryCA, size, 
                                           &evolved_steps);

    // free the 1D elementary CAs within the 2D array individually
    free_elementaryCA(all_elementaryCA, &evolved_steps); 

    /* This following 2 lines code was taken from slide 19 from the 
    lec07.pdf file in Week 8, accessed via LMS on 10/10/2024 */

    // free the memory allocated to the 2D array 
    free(all_elementaryCA);
    all_elementaryCA = NULL;

    return EXIT_SUCCESS;        // algorithms are fun!!!
}

/* In stage 0, we determine the size of the elementary cellular automata (CA)
 and calculate the updated states for the given rule and store the input 
 elementary CA in an array */
cells_t* print_stage0_output(neighbourhood_update_t* neighbourhoods, int* size){

    printf(SDELIM, 0); // print the stage delimiter for stage 0 output

    int rule; 

    scanf("%d %d", size, &rule); // recall size is already a pointer
    my_getchar(); // flush the pending new line character in input on line 2

    printf("SIZE: %d\n", *size);
    printf("RULE: %d\n", rule);

    printf(MDELIM);

    // Determine updated states for each neighbourhood using the selected rule
    determine_update_states(neighbourhoods, rule);
    print_rule(neighbourhoods);

    printf(MDELIM);
    
    cells_t * elementaryCA = malloc_elementaryCA(*size);

    // We store the input elementary CA into the array just created
    read_elementaryCA(elementaryCA, *size);
    assert(elementaryCA);

    return elementaryCA; // return the dynamically created array for further use
}

/* Calculate the update states for all neighbourhoods for a given rule */
void determine_update_states(neighbourhood_update_t * neighbourhoods, int rule){
    
    // Firstly, we store the decimal representation for each neighbourhood
    for (int i = FIRST_INDEX; i<NBRHDS; i++) {
        convert_decimal_to_binary(&neighbourhoods[i], i);
    }

    int bit_value;

    /* To calculate the rule, we start from the last index (8-1 = 7) and
    determine whether the rule is larger than 2 raised to the power of the 
    decimal value (bit value) of the neighbourhood. If so, we decrementing the 
    rule by this bit value and set the update state to ON, otherwise OFF */
    for (int i=NBRHDS + LAST_INDEX; i >= FIRST_INDEX; i--) {
        if (rule >= (bit_value = pow(BINARY_BASE,i))) {
            neighbourhoods[i].updated_state = ON;
            rule -= bit_value;
        }
        else {
            neighbourhoods[i].updated_state = OFF;
        }
    }
    
}

/* Printing the rule, the binary representation and the update states for each
 neighbourhood in the required format */
void print_rule (neighbourhood_update_t * neighbourhoods) {

    // Printing binary representation of the neighbourhood, separated by a space
    for (int i = FIRST_INDEX; i < NBRHDS; i++){

        printf(SPACE_DELIMITER);
        for (int j = FIRST_INDEX; j < NBRHD_ELEMENTS; j++) {
            printf("%d", neighbourhoods[i].binary[j]);
        }
    }

    printf("\n  "); // the first updated states has 2 leading spaces
    
    int i;
    
    // Printing update states for each neighbourhood, separated by 3 spaces
    for (i = FIRST_INDEX; i < (NBRHDS + LAST_INDEX); i++) {
        printf("%d   ", neighbourhoods[i].updated_state);
    }
    printf("%d \n", neighbourhoods[i].updated_state);

}

/* Create an array to store elementary CA, with the first and last elements
of the elementary CA also stored at the last and first index of the array
respectively for wrapping as required in stage 1 */

cells_t * malloc_elementaryCA (int size) {

    /* This following block of code was adapted from slide 19 from the lec07.pdf
    file in Week 8, accessed via LMS on 04/10/2024, altered to store cells_t */
    cells_t * elementaryCA;
    elementaryCA = (cells_t *) malloc((size + WRAPPING) * sizeof(cells_t));
    assert(elementaryCA);

    return elementaryCA;
}

/* We can read the elementary CA provided as input with given size */
void read_elementaryCA(cells_t * elementaryCA, int size) {

    // we arrange the elementary CA from indices 1 through to size (inclusive)
    int nchars = FIRST_STATE_INDEX;

    int ch;

    /* Input each state of the elementary CA and store as binary bit (0 or 1), 
    at the correct position (from index 1 to size inclusive) */
    while ((ch = my_getchar()) != NEWLINE_CHAR && nchars <= size) {
        elementaryCA[nchars] = (ch == ONSTATE);
        nchars++;
    }
    
    wrap_array_ends(elementaryCA, size);

    // print the initial elementary CA (at time t=0)
    print_elementaryCA(elementaryCA, size, FIRST_INDEX); 

}

// An improved version of getchar(); skips carriage return characters.
// NB: Adapted version of the mygetchar() function by Alistair Moffat
int my_getchar() {
    int c;
    while ((c=getchar())==CRTRNC);          // skip carriage return characters
    return c;
}

/* Print the elementary CA at a given time, time_step*/
void print_elementaryCA(cells_t * elementaryCA, int size, int time_step) {

    printf("%4d: ", time_step); // as specified in assignment specification

    // Printing character representation for each binary bit
    for (int i = FIRST_STATE_INDEX; i <= size; i++) {
        printf("%c", convert_int_to_char(elementaryCA[i]));
    }
    
    printf("\n");
};

/* Convert the binary bits (0 or 1) to its representation as characters */
char convert_int_to_char (int state) {
    if (state == ON) {
        return ONSTATE;
    } else {
        return OFFSTATE;
    }
}

/* We store the first and last state at indexes size+1 and 0 respectively 
    for wrapping so neighbourhoods of the edge states are defined */
void wrap_array_ends(cells_t * elementaryCA, int size) {

    elementaryCA[FIRST_INDEX] = elementaryCA[size];
    elementaryCA[size + NEXT_INDEX] = elementaryCA[FIRST_STATE_INDEX];
}

/* Convert the decimal representation of the neighbourhoods (in the range from 
0 to 7) to a binary representation of 3 bits */
void convert_decimal_to_binary(neighbourhood_update_t * neighbourhood, 
                               int decimal) {

    /* Again starting from last index (2) and checking if decimal is larger than
    bit value, storing correct bit, and decrementing by bit value (if needed)*/
    for (int i= NBRHD_ELEMENTS + LAST_INDEX; i >= FIRST_INDEX; i--) {

        // Store binary with the largest bit at the first index.
        if (decimal >= pow(BINARY_BASE,i)) {
            neighbourhood -> binary[NBRHD_ELEMENTS + LAST_INDEX - i] = ON;
            decimal -= pow(BINARY_BASE,i);
        }

        else {
            neighbourhood -> binary[NBRHD_ELEMENTS + LAST_INDEX - i] = OFF;   
        }
    }
}

/* In stage 1, we create a 2D array of elementary CAs (pointers) that stores how
 elementary CAs evolve over time from t=0 to t=time_step. */

cells_t ** print_stage1_output(neighbourhood_update_t * neighbourhoods, 
                        cells_t * elementaryCA, int size, int * time_step) {

    printf(SDELIM, 1); // stage delimiter for stage 1 output
    
    scanf("%d", time_step);

    /* The 2D structure will include elementary CAs at time t=0 (hence +1)
    and future elementary CAs */
    cells_t ** all_elementaryCA;
    int total_rows = *time_step + HEADER_ROW;

    /* This following 2 lines code was adapted from slide 19 from the lec07.pdf
    file in Week 8, accessed via LMS on 10/10/2024, altered to store cells_t* */
    all_elementaryCA = (cells_t **) malloc (total_rows * sizeof(cells_t *));
    assert(all_elementaryCA);

    // assign individual elementary CAs to the correct position in the 2D array
    all_elementaryCA[FIRST_INDEX] = elementaryCA;
    add_elementaryCA_to_array(all_elementaryCA, *time_step, neighbourhoods,  
                              size, FIRST_INDEX);

    // For stage 1, print number of on and off cells with required formatting
    print_on_off_state_counts(all_elementaryCA, *time_step);

    return all_elementaryCA;
}

/* We iteratively add the evolved elementary CAs to the 2D array */
void add_elementaryCA_to_array(cells_t ** all_elementaryCA, int time_step, 
        neighbourhood_update_t * neighbourhoods, int size, int starting_idx) {
   
    print_elementaryCA(all_elementaryCA[FIRST_INDEX], size, starting_idx);

    /* evolve the previous elementary CA (at time t = i-1) to time t=i and print
    this new evolved elementary CA (at time t=i) */
    for (int i = SECOND_INDEX; i <= time_step; i++) {

        all_elementaryCA[i] = compute_next_time_CA(neighbourhoods, 
                                all_elementaryCA[i + LAST_INDEX], size);
        print_elementaryCA(all_elementaryCA[i], size, starting_idx + i);
    }

    printf(MDELIM);
}

/* Calculates evolved elementary CA by comparing neighbourhoods of 3 elements,
 returning the new elementary CA */
cells_t * compute_next_time_CA(neighbourhood_update_t * neighbourhoods, 
                       cells_t * elementaryCA, int size) {
    
    cells_t * next_elementaryCA = malloc_elementaryCA(size);

    /* The previous, current and the next index are placed into this array to
    match neighbourhoods and determine the updated state */
    cells_t neighbourhood_cmp[NBRHD_ELEMENTS];

    // Elementary CA is stored from index 1 to size (0 & size + 1 is wrapping)
    for (int i = SECOND_INDEX; i <= size ; i++) {
        neighbourhood_cmp[FIRST_INDEX] = elementaryCA[i+PREVIOUS_INDEX];
        neighbourhood_cmp[SECOND_INDEX] = elementaryCA[i]; // current index
        neighbourhood_cmp[THIRD_INDEX] = elementaryCA[i+NEXT_INDEX];
        
        // add the updated state to the next elementary CA
        next_elementaryCA[i] = match_neighbourhoods(neighbourhoods, 
                                                       neighbourhood_cmp);
    }
    
    wrap_array_ends(next_elementaryCA, size);

    return next_elementaryCA;
}

/* Compares the provided neighbourhood (binary) when evolving the elementary CA
 to each neighbourhood, until a match is found, returning the updated state */
int match_neighbourhoods(neighbourhood_update_t * neighbourhoods, 
                            cells_t * neighbourhood_cmp) {

    for (int j = FIRST_INDEX; j<NBRHDS; j++) {

        // compare individual bits, and if same, return the updated state
        if (compare_binary_bits(neighbourhood_cmp, neighbourhoods[j].binary)) {
            return (neighbourhoods[j].updated_state);
        }
    }

    exit(EXIT_FAILURE);  //provided neighbourhood could not be matched
}

/* Compares each individual bit in binary format of two neighbourhoods */
int compare_binary_bits(cells_t * neighbourhood1, cells_t * neighbourhood2) {
    
    for (int i=FIRST_INDEX; i<NBRHD_ELEMENTS; i++) {
        if (neighbourhood1[i] != neighbourhood2[i]) {
            return DIFFERENT; 
        }
    }

    return SAME; // if each bit is same, the binary representation must be same
}

/* Traverses through the array from a specified starting index to print the 
 number of on and off states for a given cell */
void print_on_off_state_counts(cells_t ** all_elementaryCA, int time_steps) {

    int cell, start_position;
    scanf("%d,%d", &cell, &start_position); // in line with input format
    
    int on_states = INITIAL;
    for (int i = start_position; i <= time_steps; i++) {
        // The sum of the binary values is the sum of on states (since off is 0)
        on_states += all_elementaryCA[i][cell + NEXT_INDEX]; 
    }
    
    /* Off states is # rows in array (num_rows - start + 1) - # on_states 
    NB: the header in this case is the start position specified */
    int off_states = time_steps - start_position - on_states + HEADER_ROW;

    // Print the information in the required format
    printf(ON_OFF_STATE_COUNTS, on_states, off_states, cell, start_position);

};

/* In stage 2, we must solve the density classification problem (DCP) using 
 rules 184 and 232 a fixed number of times and printing the number of on and 
 off for a given cell from a given starting row */

cells_t ** print_stage2_output(cells_t ** all_elementaryCA, int size, 
                              steps_t * evolved_steps) {

    printf(SDELIM, 2); // print the stage delimiter for stage 2 output
    
    /* Number of times elementary CA needs to be evolved using rules 184 and
    232 respectively to solve DCP as described in assignment specification */
    int n = (size-2)/2; 
    int m = (size-1)/2; 

    /* we can now store the number of time steps the elementary CA will be 
    evolved using the specified rules so that we can free memory for all CAs */
    evolved_steps -> n_steps_184 = n; 
    evolved_steps -> m_steps_232 = m;

    int time_step = evolved_steps -> time_step;

    /* the total rows is given by all evolved steps (time_step + m + n) + 1 
    (for the header row also stored in 2D array) */
    int total_rows = time_step + HEADER_ROW + m + n; 
    
    /* This following 3 lines code was adapted from slide 19 from the lec07.pdf
    file in Week 8, accessed via LMS on 04/10/2024, altered to store cells_t */

    // Allocating additional space for to-be-evolved states to solve the DCP
    all_elementaryCA = (cells_t **) realloc (all_elementaryCA, 
                                             total_rows * sizeof (cells_t *));
    assert(all_elementaryCA);

    // Evolve the elementary CAs (rule 184 & 232) from correct row in 2D array
    solve_density_problem(all_elementaryCA, SOLVE_DCP_RULE1, 
                          n, time_step, size);
    solve_density_problem(all_elementaryCA, SOLVE_DCP_RULE2, 
                          m, time_step + n, size);

    // NB: whilst the total number of rows adds 1, the last index does not add 1
    print_on_off_state_counts(all_elementaryCA, time_step + m + n);

    printf(MDELIM);
    print_elementaryCA(all_elementaryCA[time_step], size, time_step);

    // we do not need to add the 1 since we are computing the index
    cells_t * last_row = all_elementaryCA[time_step + m + n];

    // the operator (<, > or =) - solution to the DCP
    char operator = find_operator(last_row);
    
    // Output the solution to the DCP in the required format
    printf(SOLUTION_DCP_FORMAT, time_step, operator);
    printf(THEEND);

    return all_elementaryCA;
}

/* To solve the DCP, we calculate the rule and evolve the elementary CAs using
 the update rule specified for the number of times required (given by nsteps)*/
void solve_density_problem(cells_t ** all_elementaryCA, int rule, int nsteps, 
                           int time_step, int size) {
    
    // Another neighborhood array for each rule (similar to the main one).
    neighbourhood_update_t update_rule_CA[NBRHDS]; 

    printf(RULE_STEPS_FORMAT, rule, nsteps);
    printf(MDELIM);

    determine_update_states(update_rule_CA, rule);
    add_elementaryCA_to_array(all_elementaryCA + time_step, nsteps, 
                              update_rule_CA, size, time_step);

}

/* To solve the DCP after evolving the elementary CAs, we only need to check
the first two characters (if they are different). If they are the same, we
need to only check the first character */
char find_operator(cells_t * last_row) {

    char operation;
    if (last_row[FIRST_INDEX] != last_row[SECOND_INDEX]) {
        operation = EQUIVALENT;

    } else if (last_row[FIRST_INDEX]) { 
        // the first character was on, so all must be on
        operation = GREATER;

    } else {
        // the first character was off, so all must be off
        operation = LESS;
    }

    return operation;
}

/* Repeatedly freeing dynamically allocated memory for all 1D elementary CAs
 within the 2-dimensional array of elementary CAs */
void free_elementaryCA(cells_t ** all_elementaryCA, steps_t * evolved_steps) {
    
    int total_evolved_steps = evolved_steps -> time_step + 
        evolved_steps -> n_steps_184 + evolved_steps -> m_steps_232;

    for(int i = FIRST_INDEX; i <= total_evolved_steps ; i++){

        /* This following 2 lines code was taken from slide 19 from the 
        lec07.pdf file in Week 8, accessed via LMS on 10/10/2024 */
        free(all_elementaryCA[i]);
        all_elementaryCA[i] = NULL; // these 1D arrays cannot be used further
    }

}

/* THE END -------------------------------------------------------------------*/