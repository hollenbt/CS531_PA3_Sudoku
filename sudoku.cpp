#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
using namespace std;

class Inferer {
    int board[9][9];
    int puzzle_no;
    bool MRV;
    int scheme;
    vector<vector<pair<int, int>>> assignments;

    // Search and Inference Stats
    int guess_count        = 0;
    int naked_single_apps  = 0;
    int hidden_single_apps = 0;
    int naked_pair_apps    = 0;
    int hidden_pair_apps   = 0;
    int naked_triple_apps  = 0;
    int hidden_triple_apps = 0;

    // With at most 9 elements in each set, the asymptotic time complexity advantage
    // of unordered_set vs. set is irrelevant. We can experiment with both to see
    // which has better time and space performance.

    // Cells values: if assigned, set will be emtpy
    vector<vector<set<int>>> cell;
    // Value locations: if assigned, set will be empty
    vector<vector<set<int>>> row;              // row[row #][value]:     set of col #'s
    vector<vector<set<int>>> col;              // col[col #][value]:     set of row #'s
    vector<vector<set<pair<int,int>>>> block;  // block[block #][value]: set of coords


    void increment(pair<int, int> &p, int j_low, int j_high) {
        ++p.second;
        if (p.second >= j_high) {
            p.second = j_low;
            ++p.first;
        }
    }

public:
    Inferer(int **b, int p_no, bool mrv, int s) : puzzle_no(p_no), MRV(mrv), scheme(s) {
        new_assignment_layer();
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j)
                board[i][j] = b[i][j];
    }

    int get_guess_count() { return guess_count; }

    void print() {
        int filled = 0;
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j)
                if (board[i][j] > 0)
                    ++filled;
        #ifdef DEBUG
            cout << endl;
            for (int i = 0; i < 9; ++i) {
                for (int j = 0; j < 9; ++j)
                    cout << ' ' << board[i][j];
                cout << endl;
            }
            cout << endl;
            cout << "Puzzle #: " << puzzle_no << endl;
            cout << "MRV? " << boolalpha << MRV << endl;
            cout << "Scheme: " << scheme << endl;
            cout << "Completion: " << fixed << setprecision(1) << (100 * filled) / 81.0 << '%' << endl;;
            cout << "Guess count: " << guess_count << endl;
            cout << "Naked single applications:  " << naked_single_apps << endl;
            cout << "Hidden single applications: " << hidden_single_apps << endl;
            cout << "Naked pair applications:    " << naked_pair_apps << endl;
            cout << "Hidden pair applications:   " << hidden_pair_apps << endl;
            cout << "Naked triple applications:  " << naked_triple_apps << endl;
            cout << "Hidden triple applications: " << hidden_triple_apps << endl;
            cout << endl;
        #else
            cout << puzzle_no << ',' << guess_count << ','
                 << boolalpha << MRV << ',' << scheme << ','
                 << fixed << setprecision(1) << 100 * filled / 81.0 << ','
                 << naked_single_apps << ',' << hidden_single_apps << ','
                 << naked_pair_apps << ',' << hidden_pair_apps << ','
                 << naked_triple_apps << ',' << hidden_triple_apps << endl;
        #endif
    }

    void new_assignment_layer() {
        assignments.push_back(vector<pair<int, int>>());
    }

    int infer() {
        switch(scheme) {
            case 3:  return up_to_triples();
            case 2:  return up_to_pairs();
            case 1:  return up_to_singles();
            default: return no_inference();
        }
    }

    set<int> get_guess_domain(pair<int, int> &p) {
        if (MRV) {
            int min_size = 10;
            for (int i = 0; i < 9; ++i)
                for (int j = 0; j < 9; ++j) {
                    int d_size = cell[i][j].size();
                    if (0 < d_size && d_size < min_size) {
                        p.first = i;
                        p.second = j;
                        min_size = d_size;
                    }
                }
            return cell[p.first][p.second];
        }
        else {
            for (int i = 0; i < 9; ++i)
                for (int j = 0; j < 9; ++j)
                    if (board[i][j] == 0) {
                        p.first = i;
                        p.second = j;
                        return cell[i][j];
                    }
        }
        return set<int>();
    }

    bool complete() {
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j)
                if (board[i][j] == 0) return false;
        return true;
    }
    
    bool verify() {
        return initialize();
    }

    void guess(pair<int, int> &p, int v) {
        #ifdef DEBUG
            cout << "Guess " << val << " for (" << var.first << ", " << var.second << ")" << endl;
        #endif
        ++guess_count;
        board[p.first][p.second] = v;
        assignments.back().push_back(p);
    }

    // Initializes variable domains for inference. Can also be used to perform
    // verification on a completed puzzle.
    bool initialize() {
        reset();
        for (int i = 0; i < 9; ++i) {
            for (int j = 0; j < 9; ++j) {
                int v = board[i][j];
                if (v > 0) {
                    if (assign(i, j, v, false) == 2) return false;
                }
            }
        }
        return true;
    }

    void reset() {
        cell  = vector<vector<set<int>>>(9, vector<set<int>>(9, set<int>({ 1 ,2, 3, 4, 5, 6, 7, 8, 9 })));
        row   = vector<vector<set<int>>>(9, vector<set<int>>(9, set<int>({ 0, 1 ,2, 3, 4, 5, 6, 7, 8 })));
        col   = vector<vector<set<int>>>(9, vector<set<int>>(9, set<int>({ 0, 1 ,2, 3, 4, 5, 6, 7, 8 })));
        block = vector<vector<set<pair<int,int>>>>(9, vector<set<pair<int, int>>>(9, set<pair<int, int>>()));
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j) {
                pair<int, int> p(i, j);
                int b = get_block(p);
                for (int k = 0; k < 9; ++k)
                    block[b][k].insert(p);
            }
    }

    int get_block(const pair<int, int> &p) {
        return 3 * (p.first / 3) + (p.second / 3);
    }

    /***********************
     * Variable Assignment *
     ***********************/
     // Return values:
     //   1 - successful assignment
     //   2 - conflict (i.e. empty domain for unassigned variable), need to backtrack
    int assign(int i, int j, int v, bool remember = true) {
        // Check whether assignment is valid
        if (cell[i][j].find(v) == cell[i][j].end())
            return 2;

        pair<int, int> p(i, j);
        int b = get_block(p);

        // Make assignment on board (if not during initialization)
        if (remember) {
            board[i][j] = v;
            assignments.back().push_back(p);
        }

        // Clear the sets to indicate assignment
        // (and prevent domain restriction from
        // returning empty domain error)
        cell[i][j].clear();
        row[i][v - 1].clear();
        col[j][v - 1].clear();
        block[b][v - 1].clear();

        // Remove value from rest of row and column.
        // Remove coordinates from all other values.
        for (int k = 0; k < 9; ++k) {
            if (eliminate(i, k, v)) return 2;
            if (eliminate(k, j, v)) return 2;
            if (eliminate(i, j, k + 1)) return 2;
        }

        // Remove value from rest of block
        int bi = 3 * (i / 3);
        int bj = 3 * (j / 3);
        for (int x = bi; x < bi + 3; ++x)
            for (int y = bj; y < bj + 3; ++y)
                if(eliminate(x, y, v)) return 2;

        return 1;
    }

    // Reverse all assignments made during inference
    void unassign() {
        for (const pair<int, int> &p : assignments.back())
            board[p.first][p.second] = 0;
        assignments.pop_back();
    }

    /**********************
     * Domain restriction *
     **********************/
     // Return: true if restriction makes a domain empty,
     //         false otherwise
    bool eliminate(int i, int j, int v) {
        if (cell[i][j].erase(v) && cell[i][j].empty()) return true;
        if (row[i][v - 1].erase(j) && row[i][v - 1].empty()) return true;
        if (col[j][v - 1].erase(i) && col[j][v - 1].empty()) return true;
        pair<int, int> p(i, j);
        int b = get_block(p);
        if (block[b][v - 1].erase(p) && block[b][v - 1].empty()) return true;
        return false;
    }

    /*******************
     * Inference Rules *
     *******************/
     // Return values:
     //   0 - no application found
     //   1 - application found and applied
     //   2 - conflict, need to backtrack

    int naked_single() {
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j)
                if (cell[i][j].size() == 1) {
                    #ifdef DEBUG
                        cout << "Naked single found in cell (" << i << ", " << j << ')' << endl;
                    #endif
                    ++naked_single_apps;
                    return assign(i, j, *cell[i][j].begin());
                }
        return 0;
    }

    int hidden_single() {
        for (int x = 0; x < 9; ++x) {
            for (int v = 0; v < 9; ++v) {
                if (row[x][v].size() == 1) {
                    #ifdef DEBUG
                        cout << "Hidden single found in cell (" << x << ", " << *row[x][v].begin() << ')' << endl;
                    #endif
                    ++hidden_single_apps;
                    return assign(x, *row[x][v].begin(), v + 1);
                }
                if (col[x][v].size() == 1) {
                    #ifdef DEBUG
                        cout << "Hidden single found in cell (" << *col[x][v].begin() << ", " << x << ')' << endl;
                    #endif
                    ++hidden_single_apps;
                    return assign(*col[x][v].begin(), x, v + 1);
                }
                if (block[x][v].size() == 1) {
                    pair<int, int> p = *block[x][v].begin();
                    #ifdef DEBUG
                        cout << "Hidden single found in cell (" << p.first << ", " << p.second << ')' << endl;
                    #endif
                    ++hidden_single_apps;
                    return assign(p.first, p.second, v + 1);
                }
            }
        }
        return 0;
    }

    // If we are performing this, there are no naked/hidden singles, so all row/col/block.size() > 1
    // (if referenced by a cell, otherwise could be 0). This means we only need to check that the
    // larger cell set size is > 2, and not that the smaller cell set size is >= 2.
    int naked_pair() {
        // Check by row
        for (int r = 0; r < 9; ++r) {
            for (int c1 = 0; c1 < 8; ++c1) {
                // Check for cells that have only 2 values
                if (cell[r][c1].size() == 2) {
                    set<int>::iterator it = cell[r][c1].begin();
                    bool second = false;
                    int n1 = row[r][*it - 1].size();
                    ++it;
                    int n2 = row[r][*it - 1].size();
                    if (n1 > n2) {
                        swap(n1, n2);
                        second = true;
                    }
                    if (n2 > 2) { // One of the values must appear in more than 2 cells in row (otherwise we will repeatedly find same pair)
                        it = cell[r][c1].begin();
                        if (second) ++it; // Iterate through smaller of the two column sets
                        for (set<int>::iterator cit = row[r][*it - 1].upper_bound(c1); cit != row[r][*it - 1].end(); ++cit)
                            if (cell[r][c1] == cell[r][*cit]) { // If naked pair found
                                #ifdef DEBUG
                                    cout << "Row: Naked pair found in cells (" << r << ", " << c1 << ')';
                                    cout << " and (" << r << ", " << *cit << ")" << endl;
                                #endif
                                // Eliminate the pair of values from all other cells in the row
                                for (int v : cell[r][c1]) {
                                    for (int c : row[r][v - 1])
                                        if (c != c1 && c != *cit) {
                                            #ifdef DEBUG
                                                cout << "Eliminating (" << r << ", " << c << ", " << v << ")" << endl;
                                            #endif
                                            if (eliminate(r, c, v)) return 2;
                                        }
                                }
                                ++naked_pair_apps;
                                return 1;
                            }
                    }
                }
            }
        }

        // Check by column
        for (int c = 0; c < 9; ++c) {
            for (int r1 = 0; r1 < 8; ++r1) {
                // Check for cells that have only 2 values
                if (cell[r1][c].size() == 2) {
                    set<int>::iterator it = cell[r1][c].begin();
                    bool second = false;
                    int n1 = col[c][*it - 1].size();
                    ++it;
                    int n2 = col[c][*it - 1].size();
                    if (n1 > n2) {
                        swap(n1, n2);
                        second = true;
                    }
                    if (n2 > 2) { // One of the values must appear in more than 2 cells in column (otherwise we will repeatedly find same pair)
                        it = cell[r1][c].begin();
                        if (second) ++it; // Iterate through smaller of the two row sets
                        for (set<int>::iterator rit = col[c][*it - 1].upper_bound(r1); rit != col[c][*it - 1].end(); ++rit)
                            if (cell[r1][c] == cell[*rit][c]) { // If naked pair found
                                #ifdef DEBUG
                                    cout << "Column: Naked pair found in cells (" << r1 << ", " << c << ')';
                                    cout << " and (" << *rit << ", " << c << ")" << endl;
                                #endif
                                // Eliminate the pair of values
                                for (int v : cell[r1][c]) {
                                    // From all other cells in the column
                                    for (int r : col[c][v - 1])
                                        if (r != r1 && r != *rit) {
                                            #ifdef DEBUG
                                                cout << "Eliminating (" << r << ", " << c << ", " << v << ")" << endl;
                                            #endif
                                            if (eliminate(r, c, v)) return 2;
                                        }
                                }
                                ++naked_pair_apps;
                                return 1;
                            }
                    }
                }
            }
        }

        // Check by block
        for (int b = 0; b < 9; ++b) {
            int r_beg = 3 * (b / 3);
            int c_beg = 3 * (b % 3);
            for (int r1 = r_beg; r1 < r_beg + 3; ++r1) {
                for (int c1 = c_beg; c1 < c_beg + 3; ++c1) {
                    // Check for cells that have only 2 values
                    if (cell[r1][c1].size() == 2) {
                        pair<int, int> p1(r1, c1);
                        set<int>::iterator it = cell[r1][c1].begin();
                        bool second = false;
                        int n1 = block[b][*it - 1].size();
                        ++it;
                        int n2 = block[b][*it - 1].size();
                        if (n1 > n2) {
                            swap(n1, n2);
                            second = true;
                        }
                        if (n2 > 2) { // One of the values must appear in more than 2 cells in block (otherwise we will repeatedly find same pair)
                            it = cell[r1][c1].begin();
                            if (second) ++it; // Iterate through smaller of the two coordinate sets
                            for (set<pair<int, int>>::iterator pit = block[b][*it - 1].upper_bound(p1); pit != block[b][*it - 1].end(); ++pit)
                                if (cell[r1][c1] == cell[pit->first][pit->second]) { // If naked pair found
                                    #ifdef DEBUG
                                        cout << "Block: Naked pair found in cells (" << r1 << ", " << c1 << ')';
                                        cout << " and (" << pit->first << ", " << pit->second << ")" << endl;
                                    #endif
                                    // Eliminate the pair of values from all other cells in the block
                                    for (int v : cell[r1][c1]) {
                                        for (const pair<int, int> &p : block[b][v - 1])
                                            if (p != p1 && p != *pit) {
                                                #ifdef DEBUG
                                                    cout << "Eliminating (" << p.first << ", " << p.second << ", " << v << ")" << endl;
                                                #endif
                                                if (eliminate(p.first, p.second, v)) return 2;
                                            }
                                    }
                                    ++naked_pair_apps;
                                    return 1;
                                }
                        }
                    }
                }
            }
        }

        return 0;
    }

    // If we are performing this, there are no naked/hidden singles or naked pairs, so
    // we can ignore pairs of cells where the larger value set size is 2 or less.
    int hidden_pair() {
        for (int x = 0; x < 9; ++x) {
            for (int v = 0; v < 8; ++v) {
                // Check row for values that appear in only 2 cells
                if (row[x][v].size() == 2) {
                    set<int>::iterator it = row[x][v].begin();
                    bool second = false;
                    int n1 = cell[x][*it].size();
                    ++it;
                    int n2 = cell[x][*it].size();
                    if (n1 > n2) {
                        swap(n1, n2);
                        second = true;
                    }
                    if (n2 > 2) { // One of the cells must have more than 2 values (otherwise it has already been handled as a naked pair)
                        it = row[x][v].begin();
                        if (second) ++it; // Iterate through smaller of the two value sets
                        for (set<int>::iterator vit = cell[x][*it].upper_bound(v + 1); vit != cell[x][*it].end(); ++vit)
                            if (row[x][v] == row[x][*vit - 1]) { // If hidden pair found
                                #ifdef DEBUG
                                    cout << "Row: Hidden pair found with values " << v + 1 << ", " << *vit << endl;
                                #endif
                                // Eliminate all other values from the pair of cells.
                                for (int c : row[x][v]) {
                                    for (int val : cell[x][c])
                                        if (val != v + 1 && val != *vit) {
                                            #ifdef DEBUG
                                                cout << "Eliminating (" << x << ", " << c << ", " << val << ")" << endl;
                                            #endif
                                            if (eliminate(x, c, val)) return 2;
                                        }
                                }
                                ++hidden_pair_apps;
                                return 1;
                            }
                    }
                }

                // Check column for values that appear in only 2 cells
                if (col[x][v].size() == 2) {
                    set<int>::iterator it = col[x][v].begin();
                    bool second = false;
                    int n1 = cell[*it][x].size();
                    ++it;
                    int n2 = cell[*it][x].size();
                    if (n1 > n2) {
                        swap(n1, n2);
                        second = true;
                    }
                    if (n2 > 2) { // One of the cells must have more than 2 values (otherwise it has already been handled as a naked pair)
                        it = col[x][v].begin();
                        if (second) ++it; // Iterate through smaller of the two value sets
                        for (set<int>::iterator vit = cell[*it][x].upper_bound(v + 1); vit != cell[*it][x].end(); ++vit)
                            if (col[x][v] == col[x][*vit - 1]) { // If hidden pair found
                                #ifdef DEBUG
                                    cout << "Column: Hidden pair found with values " << v + 1 << ", " << *vit << endl;
                                #endif
                                // Eliminate all other values from the pair of cells.
                                for (int r : col[x][v]) {
                                    for (int val : cell[r][x])
                                        if (val != v + 1 && val != *vit) {
                                            #ifdef DEBUG
                                                cout << "Eliminating (" << r << ", " << x << ", " << val << ")" << endl;
                                            #endif
                                            if (eliminate(r, x, val)) return 2;
                                        }
                                }
                                ++hidden_pair_apps;
                                return 1;
                            }
                    }
                }

                // Check block for values that appear in only 2 cells
                if (block[x][v].size() == 2) {
                    set<pair<int, int>>::iterator it = block[x][v].begin();
                    bool second = false;
                    int n1 = cell[it->first][it->second].size();
                    ++it;
                    int n2 = cell[it->first][it->second].size();
                    if (n1 > n2) {
                        swap(n1, n2);
                        second = true;
                    }
                    if (n2 > 2) { // One of the cells must have more than 2 values (otherwise it has already been handled as a naked pair)
                        it = block[x][v].begin();
                        if (second) ++it; // Iterate through smaller of the two value sets
                        for (set<int>::iterator vit = cell[it->first][it->second].upper_bound(v + 1); vit != cell[it->first][it->second].end(); ++vit)
                            if (block[x][v] == block[x][*vit - 1]) { // If hidden pair found
                                #ifdef DEBUG
                                    cout << "Block: Hidden pair found with values " << v + 1 << ", " << *vit << endl;
                                #endif
                                // Eliminate all other values from the pair of cells.
                                for (const pair<int, int> &p : block[x][v]) {
                                    for (int val : cell[p.first][p.second])
                                        if (val != v + 1 && val != *vit) {
                                            #ifdef DEBUG
                                                cout << "Eliminating (" << p.first << ", " << p.second << ", " << val << ")" << endl;
                                            #endif
                                            if (eliminate(p.first, p.second, val)) return 2;
                                        }
                                }
                                ++hidden_pair_apps;
                                return 1;
                            }
                    }
                }
            }
        }
        return 0;
    }

    // No naked/hidden singles, no naked/hidden pairs if we are performing this.
    // First, no cell in the triple can have more than 3 elements. Second, the
    // pairwise union of any two cells in the triple must contain exactly 3 values
    // (because of the inference rule application precedence). These facts give us
    // ways to prune the space of triples to consider.
    int naked_triple() {
        // Check by row
        for (int r = 0; r < 9; ++r) {
            for (int c1 = 0; c1 < 7; ++c1) {
                int s1 = cell[r][c1].size();
                if (0 < s1 && s1 < 4) { // note: because no naked singles, 1 is not a possible size for any cell
                    for (int c2 = c1 + 1; c2 < 8; ++c2) {
                        int s2 = cell[r][c2].size();
                        if (0 < s2 && s2 < 4) {
                            set<int> union12(cell[r][c1]);
                            union12.insert(cell[r][c2].begin(), cell[r][c2].end());
                            if (union12.size() == 3) {
                                for (int c3 = c2 + 1; c3 < 9; ++c3) {
                                    int s3 = cell[r][c3].size();
                                    if (0 < s3 && s3 < 4) {
                                        set<int> union123(union12);
                                        union123.insert(cell[r][c3].begin(), cell[r][c3].end());
                                        if (union123.size() == 3) {
                                            // At this point, before declaring it a naked triple, we must check that one of the values
                                            // appears in a cell in the row outside of the triple (otherwise, we will repeatedly find
                                            // the same triple).
                                            #ifdef DEBUG
                                                cout << "Row: Naked triple found in cells (" << r << ", " << c1 << "), (";
                                                cout << r << ", " << c2 << "), and (" << r << ", " << c3 << ")" << endl;
                                            #endif
                                            // Eliminate the triple of values from all other cells in the row
                                            bool eliminated = false;
                                            for (int v : union123) {
                                                for (int c : row[r][v - 1])
                                                    if (c != c1 && c != c2 && c != c3) {
                                                        eliminated = true;
                                                        #ifdef DEBUG
                                                            cout << "Eliminating (" << r << ", " << c << ", " << v << ")" << endl;
                                                        #endif
                                                        if (eliminate(r, c, v)) return 2;
                                                    }
                                            }
                                            if (eliminated) {
                                                ++naked_triple_apps;
                                                return 1;
                                            }
                                            #ifdef DEBUG
                                                else cout << "Nothing eliminated." << endl;
                                            #endif
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Check by column
        for (int c = 0; c < 9; ++c) {
            for (int r1 = 0; r1 < 7; ++r1) {
                int s1 = cell[r1][c].size();
                if (0 < s1 && s1 < 4) { // note: because no naked singles, 1 is not a possible size for any cell
                    for (int r2 = r1 + 1; r2 < 8; ++r2) {
                        int s2 = cell[r2][c].size();
                        if (0 < s2 && s2 < 4) {
                            set<int> union12(cell[r1][c]);
                            union12.insert(cell[r2][c].begin(), cell[r2][c].end());
                            if (union12.size() == 3) {
                                for (int r3 = r2 + 1; r3 < 9; ++r3) {
                                    int s3 = cell[r3][c].size();
                                    if (0 < s3 && s3 < 4) {
                                        set<int> union123(union12);
                                        union123.insert(cell[r3][c].begin(), cell[r3][c].end());
                                        if (union123.size() == 3) {
                                            // At this point, before declaring it a naked triple, we must check that one of the values
                                            // appears in a cell in the column outside of the triple (otherwise, we will repeatedly find
                                            // the same triple).
                                            #ifdef DEBUG
                                                cout << "Column: Naked triple found in cells (" << r1 << ", " << c << "), (";
                                                cout << r2 << ", " << c << "), and (" << r3 << ", " << c << ")" << endl;
                                            #endif
                                            // Eliminate the triple of values from all other cells in the column
                                            bool eliminated = false;
                                            for (int v : union123) {
                                                for (int r : col[c][v - 1])
                                                    if (r != r1 && r != r2 && r != r3) {
                                                        eliminated = true;
                                                        #ifdef DEBUG
                                                            cout << "Eliminating (" << r << ", " << c << ", " << v << ")" << endl;
                                                        #endif
                                                        if (eliminate(r, c, v)) return 2;
                                                    }
                                            }
                                            if (eliminated) {
                                                ++naked_triple_apps;
                                                return 1;
                                            }
                                            #ifdef DEBUG
                                                else cout << "Nothing eliminated." << endl;
                                            #endif
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Check by block
        for (int b = 0; b < 9; ++b) {
            int r_beg = 3 * (b / 3);
            int c_beg = 3 * (b % 3);
            pair<int, int> p1_end(r_beg + 2, c_beg + 1);
            pair<int, int> p2_end(r_beg + 2, c_beg + 2);
            pair<int, int> p3_end(r_beg + 3, c_beg);
            for (pair<int, int> p1(r_beg, c_beg); p1 < p1_end; increment(p1, c_beg, c_beg + 3)) {
                int s1 = cell[p1.first][p1.second].size();
                if (0 < s1 && s1 < 4) { // note: because no naked singles, 1 is not a possible size for any cell
                    pair<int, int> p2(p1);
                    for (increment(p2, c_beg, c_beg + 3); p2 < p2_end; increment(p2, c_beg, c_beg + 3)) {
                        int s2 = cell[p2.first][p2.second].size();
                        if (0 < s2 && s2 < 4) {
                            set<int> union12(cell[p1.first][p1.second]);
                            union12.insert(cell[p2.first][p2.second].begin(), cell[p2.first][p2.second].end());
                            if (union12.size() == 3) {
                                pair<int, int> p3(p2);
                                for (increment(p3, c_beg, c_beg + 3); p3 < p3_end; increment(p3, c_beg, c_beg + 3)) {
                                    int s3 = cell[p3.first][p3.second].size();
                                    if (0 < s3 && s3 < 4) {
                                        set<int> union123(union12);
                                        union123.insert(cell[p3.first][p3.second].begin(), cell[p3.first][p3.second].end());
                                        if (union123.size() == 3) {
                                            // At this point, before declaring it a naked triple, we must check that one of the values
                                            // appears in a cell in the block outside of the triple (otherwise, we will repeatedly find
                                            // the same triple).
                                            #ifdef DEBUG
                                                cout << "Block: Naked triple found in cells (" << p1.first << ", " << p1.second << "), (";
                                                cout << p2.first << ", " << p2.second << "), and (" << p3.first << ", " << p3.second << ")" << endl;
                                            #endif
                                            // Eliminate the triple of values from all other cells in the block
                                            bool eliminated = false;
                                            for (int v : union123) {
                                                for (const pair<int, int> &p : block[b][v - 1])
                                                    if (p != p1 && p != p2 && p != p3) {
                                                        eliminated = true;
                                                        #ifdef DEBUG
                                                            cout << "Eliminating (" << p.first << ", " << p.second << ", " << v << ")" << endl;
                                                        #endif
                                                        if (eliminate(p.first, p.second, v)) return 2;
                                                    }
                                            }
                                            if (eliminated) {
                                                ++naked_triple_apps;
                                                return 1;
                                            }
                                            #ifdef DEBUG
                                                else cout << "Nothing eliminated." << endl;
                                            #endif
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return 0;
    }

    // If we are performing this, there are no naked/hidden singles or pairs, and
    // no naked triples. First, no value in the triple can appear in more than 3 cells.
    // Second, the pairwise union of any two values in the triple must contain exactly
    // 3 cells (because of the inference rule application precedence). These facts give
    // us ways to prune the space of triples to consider.
    int hidden_triple() {
        for (int x = 0; x < 9; ++x) {
            for (int v1 = 0; v1 < 7; ++v1) {
                // Check row for values that appear in at most 3 cells
                int s1 = row[x][v1].size();
                if (0 < s1 && s1 < 4) { // note: because no hidden singles, 1 is not a possible size for any value
                    for (int v2 = v1 + 1; v2 < 8; ++v2) {
                        int s2 = row[x][v2].size();
                        if (0 < s2 && s2 < 4) {
                            set<int> union12(row[x][v1]);
                            union12.insert(row[x][v2].begin(), row[x][v2].end());
                            if (union12.size() == 3) {
                                for (int v3 = v2 + 1; v3 < 9; ++v3) {
                                    int s3 = row[x][v3].size();
                                    if (0 < s3 && s3 < 4) {
                                        set<int> union123(union12);
                                        union123.insert(row[x][v3].begin(), row[x][v3].end());
                                        if (union123.size() == 3) {
                                            // At this point, before declaring it a hidden triple, we must check that one of the cells
                                            // in the triple has some value outside the triple (otherwise, it has already been handled
                                            // as a naked triple).
                                            #ifdef DEBUG
                                                cout << "Row: Hidden triple found with values " << v1 + 1 << ", " << v2 + 1 << ", and " << v3 + 1 << endl;
                                            #endif
                                            // Eliminate all other values from the triple of cells.
                                            bool eliminated = false;
                                            for (int c : union123) {
                                                for (int v : cell[x][c])
                                                    if (v != v1 + 1 && v != v2 + 1 && v != v3 + 1) {
                                                        eliminated = true;
                                                        #ifdef DEBUG
                                                            cout << "Eliminating (" << x << ", " << c << ", " << v << ")" << endl;
                                                        #endif
                                                        if (eliminate(x, c, v)) return 2;
                                                    }
                                            }
                                            if (eliminated) {
                                                ++hidden_triple_apps;
                                                return 1;
                                            }
                                            #ifdef DEBUG
                                                else cout << "Nothing eliminated." << endl;
                                            #endif
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Check column for values that appear in at most 3 cells
                s1 = col[x][v1].size();
                if (0 < s1 && s1 < 4) { // note: because no hidden singles, 1 is not a possible size for any value
                    for (int v2 = v1 + 1; v2 < 8; ++v2) {
                        int s2 = col[x][v2].size();
                        if (0 < s2 && s2 < 4) {
                            set<int> union12(col[x][v1]);
                            union12.insert(col[x][v2].begin(), col[x][v2].end());
                            if (union12.size() == 3) {
                                for (int v3 = v2 + 1; v3 < 9; ++v3) {
                                    int s3 = col[x][v3].size();
                                    if (0 < s3 && s3 < 4) {
                                        set<int> union123(union12);
                                        union123.insert(col[x][v3].begin(), col[x][v3].end());
                                        if (union123.size() == 3) {
                                            // At this point, before declaring it a hidden triple, we must check that one of the cells
                                            // in the triple has some value outside the triple (otherwise, it has already been handled
                                            // as a naked triple).
                                            #ifdef DEBUG
                                                cout << "Column: Hidden triple found with values " << v1 + 1 << ", " << v2 + 1 << ", and " << v3 + 1 << endl;
                                            #endif
                                            // Eliminate all other values from the triple of cells.
                                            bool eliminated = false;
                                            for (int r : union123) {
                                                for (int v : cell[r][x])
                                                    if (v != v1 + 1 && v != v2 + 1 && v != v3 + 1) {
                                                        eliminated = true;
                                                        #ifdef DEBUG
                                                            cout << "Eliminating (" << r << ", " << x << ", " << v << ")" << endl;
                                                        #endif
                                                        if (eliminate(r, x, v)) return 2;
                                                    }
                                            }
                                            if (eliminated) {
                                                ++hidden_triple_apps;
                                                return 1;
                                            }
                                            #ifdef DEBUG
                                                else cout << "Nothing eliminated." << endl;
                                            #endif
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Check block for values that appear in at most 3 cells
                s1 = block[x][v1].size();
                if (0 < s1 && s1 < 4) { // note: because no hidden singles, 1 is not a possible size for any value
                    for (int v2 = v1 + 1; v2 < 8; ++v2) {
                        int s2 = block[x][v2].size();
                        if (0 < s2 && s2 < 4) {
                            set<pair<int, int>> union12(block[x][v1]);
                            union12.insert(block[x][v2].begin(), block[x][v2].end());
                            if (union12.size() == 3) {
                                for (int v3 = v2 + 1; v3 < 9; ++v3) {
                                    int s3 = block[x][v3].size();
                                    if (0 < s3 && s3 < 4) {
                                        set<pair<int, int>> union123(union12);
                                        union123.insert(block[x][v3].begin(), block[x][v3].end());
                                        if (union123.size() == 3) {
                                            // At this point, before declaring it a hidden triple, we must check that one of the cells
                                            // in the triple has some value outside the triple (otherwise, it has already been handled
                                            // as a naked triple).
                                            #ifdef DEBUG
                                                cout << "Block: Hidden triple found with values " << v1 + 1 << ", and " << v2 + 1 << ", " << v3 + 1 << endl;
                                            #endif
                                            // Eliminate all other values from the triple of cells.
                                            bool eliminated = false;
                                            for (const pair<int, int> &p : union123) {
                                                for (int v : cell[p.first][p.second])
                                                    if (v != v1 + 1 && v != v2 + 1 && v != v3 + 1) {
                                                        eliminated = true;
                                                        #ifdef DEBUG
                                                            cout << "Eliminating (" << p.first << ", " << p.second << ", " << v << ")" << endl;
                                                        #endif
                                                        if (eliminate(p.first, p.second, v)) return 2;
                                                    }
                                            }
                                            if (eliminated) {
                                                ++hidden_triple_apps;
                                                return 1;
                                            }
                                            #ifdef DEBUG
                                                else cout << "Nothing eliminated." << endl;
                                            #endif
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return 0;
    }

    /*********************
     * Inference Schemes *
     *********************/
    int no_inference() {
        return 0;
    }

    int up_to_singles() {
        int r = 1;
        while (r == 1) {
            r = naked_single();
            if (r == 0)
                r = hidden_single();
        }
        return r;
    }

    int up_to_pairs() {
        int r = 1;
        while (r == 1) {
            r = naked_single();
            if (r == 0) {
                r = hidden_single();
                if (r == 0) {
                    r = naked_pair();
                    if (r == 0)
                        r = hidden_pair();
                }
            }
        }
        return r;
    }

    int up_to_triples() {
        int r = 1;
        while (r == 1) {
            r = naked_single();
            if (r == 0) {
                r = hidden_single();
                if (r == 0) {
                    r = naked_pair();
                    if (r == 0) {
                        r = hidden_pair();
                        if (r == 0) {
                            r = naked_triple();
                            if (r == 0)
                                r == hidden_triple();
                        }
                    }
                }
            }
        }
        return r;
    }
};

/***********************
 * Backtracking Search *
 ***********************/
bool backtracking_search(Inferer &inferer) {
    #ifdef DEBUG
        cout << "New stack frame" << endl;
    #endif
    // If the last guess created a conflict or if a conflict
    // is encountered during inference, backtrack immediately
    if (!inferer.initialize() || inferer.infer() == 2) {
        #ifdef DEBUG
            cout << "Backtracking..." << endl;
        #endif
        return false;
    }
    // If puzzle is solved, hooray!
    if (inferer.complete()) return true;

    pair<int, int> var;
    set<int> domain = inferer.get_guess_domain(var);
    for (int val : domain) {
        inferer.new_assignment_layer();
        if (inferer.get_guess_count() >= 1000) return true;
        inferer.guess(var, val);
        if (backtracking_search(inferer)) return true;
        inferer.unassign();
    }
    // If no values in the domain lead to a solution, we need to backtrack further
    #ifdef DEBUG
        cout << "Backtracking..." << endl;
    #endif
    return false;
}

void read_puzzle(vector<int**> &board) {
    board.push_back(new int*[9]);
    int **b = board.back();

    for (int i = 0; i < 9; ++i)
        b[i] = new int[9];

    string s;
    getline(cin, s); // difficulty level is read in this line, if needed for something
    for (int i = 0; i < 9; ++i) {
        getline(cin, s);
        b[i][0] =  s[0] - '0';
        b[i][1] =  s[1] - '0';
        b[i][2] =  s[2] - '0';
        b[i][3] =  s[4] - '0';
        b[i][4] =  s[5] - '0';
        b[i][5] =  s[6] - '0';
        b[i][6] =  s[8] - '0';
        b[i][7] =  s[9] - '0';
        b[i][8] = s[10] - '0';
    }
    getline(cin, s);
}

int main() {
    vector<int**> board;

    while (board.size() < 77)
        read_puzzle(board);

    for (bool m : { true, false }) {
        for (int s : { 0, 1, 2, 3 }) { 
            for (int i = 0; i < board.size(); ++i) {
                Inferer inferer(board[i], i + 1, m, s);
                backtracking_search(inferer);
                inferer.print();
            }
        }
    }

    for (int **p : board) {
        for (int i = 0; i < 9; ++i)
            delete[] p[i];
        delete[] p;
    }

    return 0;
}