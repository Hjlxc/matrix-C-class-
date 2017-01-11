#ifndef MY_MATRIX_H_
#define MY_MATRIX_H_

#include<iostream>
#include<vector>
#include<utility>
#include<iomanip>
#include<complex>
#include <stdlib.h>     
using std::cout;
using std::endl;
using std::vector;
using std::complex;

typedef unsigned int uint;

// Forward declarations.
template<typename entry_element> class matrix;

//non-mumber function

//generate a matrix of size row_num*col_num
//all the entry are satisfy Uniform distribution between provided begin(default=0) and end(default=1)
//the entry are generalized randomly using function rand()
template<typename entry_element>
matrix<entry_element> urand(int row_num,int col_num=1,const entry_element& begin = entry_element(), const entry_element& end = entry_element(1));
//reshape a vector to a matrix with given row_num and col_num
//if the size of vector is not equal to the size of the reshaped matrix, an error message will printed
template<typename entry_element>
matrix<entry_element> reshape(const vector<entry_element>& vect, int row_num, int col_num);

//class defination
template <typename entry_element> 
class matrix {
private:
	vector<vector<entry_element>> my_matrix;
	uint rownum; //num_of_row
	uint colnum; //num_of_col
	//private mumber function
	void col_copy(uint target, uint source); //copy the comtain of source'th col to target's col
	//diagonalization the matrix using Gaussian elemination, used to find the inverse of a matrix
	//the colnum should larger than rownum
	//return true if the matrix can be diagonalized, otherwise false
	bool Gaussian_Diagonalization();
public:

	matrix() : rownum(0), colnum(0),my_matrix() {}  //default construct

    //construct a matrix with a 2-D vector, if the nested vector has different size
	//zeros will be added at the end
	matrix(const vector<vector<entry_element>> & m);  //constructor 1

	//construct a vector(1-D matrix) with a 1-D vector, the matrix will have the same row as v and only 1 col
	matrix(const vector<entry_element> & v); //constructor 2

    //construct a matrix with 'row' rows and 'col' columns and fill with entry,
	//row and col must > 0, otherwise call the default construct
	matrix(const int row,const int col, const entry_element & entry); //constructor 3

	//construct a diagonial matrix with size 'n*n' and fill diagonal entry with entry,
	//n must >0, otherwise call the default construct
	matrix(const int n, const entry_element & entry); //constructor 4

	//explicit type convert
	template<typename T>
	explicit matrix(const matrix<T> & m);

	//used auto-generated copy construct
	//matrix(const matrix<entry_element> & m);

	//basic number function
	uint row_num() const { return rownum; }
	uint col_num() const { return colnum; }
	void show() const; //display the matrix
	bool is_empty() const { return rownum==0||colnum==0; } //a matrix without any entry is an empty matrix
	void clear() { my_matrix.clear();rownum = colnum = 0; }
	vector<entry_element> matrix_to_vector() const; //convert matrix to a vector, the elements will be order by row
	//resize matrix with the new rownum and colnum number
	//if new rownum <= rownum, then truncate matrix with the new row number
	//else add extra rows until meet the new rownum, the added rows will filled with given value(default=0)
	//same for the column
	void resize(int new_rownum, int new_colnum, const entry_element & value=entry_element());
	/////////////////////////////////////////////////////
	//operator overload
	/////////////////////////////////////////////////////
	//addition of two matrix with same size *this+m
	template<typename T>
	matrix<entry_element> operator+(const matrix<T>& m) const; //matrix must have same size, otherwise show a warning and return a enpty matrix
	//substitution two matrix with same size *this-m
	template<typename T>
	matrix<entry_element> operator-(const matrix<T>& m) const; //matrix must have same size, otherwise show a warning and return a enpty matrix
	//multiplication two matrixx (*this)*m
	template<typename T>
	matrix<entry_element> operator*(const matrix<T>& m) const; //this->colnum == m.rownum, otherwise show a warning and return a enpty matrix
	//multiplication matrix with scalar
	template<typename T>
	matrix<entry_element> operator*(const T& scalar) const; 
	//matrix devided my a scalar
	template<typename T>
	matrix<entry_element> operator/(const T& scalar) const;
	//sovle the linear system Ax=b using Gaussian elemination, where A must be square and x and b are vectors
	matrix<double> operator%(const matrix<entry_element> & m)const;
	// negative -(*this)
	matrix<entry_element> operator-() const; 
	// transpose
	matrix<entry_element> operator~() const;
	// access the entry
	const vector<entry_element>& operator[](uint row) const { return my_matrix[row]; }
	// access and change the entry
	vector<entry_element>& operator[](uint row) { return my_matrix[row]; }
	
	
	/////////////////////////////////////////////////////
	//matrix operator
	/////////////////////////////////////////////////////

	// add an addition row at back with all elements fill with given value(default=0)
	void add_row(const entry_element& value = entry_element()); //type1

	// add an addition row at a given location, the added row are filled with given value
	// if location<0 || location>rownum, print a error message and do nothing
	// else add at corresponding location (location==0 add at begining, location==rownum add at end)
	void add_row(int location, const entry_element& value); //type2

	//add an addition row at back with a given vector
	//if size of the v is larger than colnum, the v will be truncated to the correct size
	//if size of the v is smaller than colnum, 0 will be added to the left entry
	void add_row(const vector<entry_element>& v);//type3

	//add an addition row at given location with a given vector
	//the affect of the location and v are the same as type2 and type3
	void add_row(int location, const vector<entry_element>& v);//type4

	//remove a row of the given location
	// if location<0 ||location>rownum-1 print a error message and do nothing
	// else remove the corresponding location
	void remove_row(int location);


	// similiar as add_row type1
	void add_col(const entry_element& value = entry_element()); //type1
	// similiar as add_row type2
	void add_col(int location, const entry_element& value); //type2
	// similiar as add_row type3
	void add_col(const vector<entry_element>& v);//type3
	// similiar as add_row type4
	void add_col(int location, const vector<entry_element>& v);//type4
	// similiar as remove_row type1
	void remove_col(int location);

	//add entry
	//row_index<0 or col_index<0 , an error will be printed
	//if row_index and col_index are all insize the matrix, the function will reset the corresponding entry to the value(default=0)
	//if row_index and/or col_index are out of matrix size, the extra row and/or col will be added until the row_index and col_index are insize matrix
	//the corresponding entry will be set to the value and rest will be zero
	void add_entry(int row_index, int col_index, const entry_element& value = entry_element());

	matrix<entry_element> get_row(int row_index) const;
	matrix<entry_element> get_col(int col_index) const;
	matrix<entry_element> submatrix(int row_index, int col_index, int num_of_row, int num_of_col);

	//friend function
	//calculate the inverse of the matrix m using gaussian elimination
	//m must be a squrse matrix, otherwise an error message will be printed and return a empty matrix
	//if m is singular, a message will be printed and return an empty matrix
	template<typename entry_element>
	friend matrix<double> inv(const matrix<entry_element> & m); //calculate the inverse of the matrix m

	template<typename T,typename entry_element>
	friend matrix<entry_element> operator*(T scalar, const matrix<entry_element>& m);
}; //class matrix
//private mumber function
template<typename entry_element>
void matrix<entry_element>::col_copy(uint target, uint source)
{
	for (uint i = 0;i < rownum;++i)
		my_matrix[i][target] = my_matrix[i][source];
}

template<typename entry_element>
bool matrix<entry_element>::Gaussian_Diagonalization()
{
	for (uint i = 0;i < rownum;++i) { //elemination step1, change the matrix to upper triangular matrix using Gaussian elemination
		if (my_matrix[i][i] == 0) { //if the pivot is zero
			for (uint j = i + 1;j < rownum;++j)  //switch it with another row which is non - zero, if there are any
				if (my_matrix[j][i] != 0) {
					swap(my_matrix[i], my_matrix[j]);
					break;
				}
		}
		if (my_matrix[i][i] == 0) //if there is no non-zero pivot
			return false;
		entry_element pivot = my_matrix[i][i];
		for (uint j = i;j <colnum;++j) //set the pivot to one and adjust all rest elements in the same row
			my_matrix[i][j] /= pivot;
		for (uint j = i + 1;j < rownum;++j) {//applied the elemination to the rest row
			if (my_matrix[j][i] == 0)
				continue;
			entry_element parameter = my_matrix[j][i] / my_matrix[i][i];
			for (uint k = i;k < colnum;++k)
				my_matrix[j][k] -= my_matrix[i][k] * parameter;
		}
	} //end of elemination step1
	//the matrix shoud become a Upper triangular matrix now
	for (uint i = 1;i < rownum;++i) { //elemination step2, diagonalization
		for (uint j = 0;j < i;++j) {
			if (my_matrix[j][i] == entry_element())
				continue;
			entry_element pivot = my_matrix[j][i];
			for (uint k = i;k < colnum;++k)
				my_matrix[j][k] -= pivot*my_matrix[i][k];
		}
	}
	return true;
}


//constructor
template<typename entry_element>
matrix<entry_element>::matrix(const vector<vector<entry_element>>& m):my_matrix(m)
{
	if (m.empty()) 
		rownum = colnum = 0;
	else {
		rownum = m.size(); //set num of row
		colnum = 0; //set num of col to zero 
		for (auto row : m)
			colnum = colnum > row.size() ? colnum : row.size(); //find the bigest col value
		for (uint i = 0;i < rownum;++i) {
			uint j = my_matrix[i].size();
			while (j++ < colnum)
				my_matrix[i].push_back(entry_element()); //add the rest entry to zero
		}
	}
}
template<typename entry_element>
matrix<entry_element>::matrix(const vector<entry_element>& v)
{
	if(v.empty())
		matrix<entry_element>();
	else {
		rownum = v.size();
		colnum = 1;
		my_matrix = vector<vector<entry_element>>(rownum, vector<entry_element>(1, entry_element()));
		for (uint i = 0;i < rownum;++i)
			my_matrix[i][0] = v[i];
	}
}
template<typename entry_element>
matrix<entry_element>::matrix(const int row,const int col,const entry_element& entry)
{
	if (row < 1 || col < 1)
		matrix<entry_element>();
	else {
		rownum = uint(row);
		colnum = uint(col);
		my_matrix = vector<vector<entry_element>>(row, vector<entry_element>(col, entry));
	}
}

template<typename entry_element>
matrix<entry_element>::matrix(const int n, const entry_element & entry)
{
	if (n < 1)
		matrix<entry_element>();
	else {
		rownum = colnum = unsigned int(n);
		my_matrix = vector<vector<entry_element>>(n, vector<entry_element>(n, entry_element()));
		for (int i = 0;i < n;++i)
			my_matrix[i][i] = entry;
	}
}

template<typename entry_element>
template<typename T>
matrix<entry_element>::matrix(const matrix<T>& m)
{
	rownum = m.row_num();
	colnum = m.col_num();
	my_matrix = vector<vector<entry_element>>(rownum, vector<entry_element>(colnum, entry_element()));
	for (uint i = 0;i < rownum;++i)
		for (uint j = 0;j < colnum;++j)
			my_matrix[i][j] = m[i][j];
}


//basic mumber functions

template<typename entry_element>
void matrix<entry_element>::show() const
{
	if (rownum*colnum == 0) {
		cout << "matrix is empty" << endl;
		return;
	}
	for (auto row : my_matrix) {
		for (auto entry : row)
			cout<<fixed << setw(8) << setprecision(2) << entry;
		cout << endl;
	}
}
template<typename entry_element>
void matrix<entry_element>::resize(int new_rownum, int new_colnum, const entry_element & value)
{
	if (new_rownum < 0 || new_colnum < 0) {
		cout << "resize error: can't resize to negitave size\n";
		return;
	}
	if (new_rownum == 0 || new_colnum == 0) {//well...
		(*this).clear(); return;
	}
	if ((*this).is_empty())
		(*this) = matrix<entry_element>(new_rownum, new_colnum, value);
	//adjust row
	if (new_rownum < rownum) {
		my_matrix.resize(new_rownum);
		rownum = new_rownum;
	}
	else if (new_rownum > rownum){
		while (rownum < new_rownum)
			(*this).add_row(value);
	}
	//adjust colnum
	if (new_colnum < colnum) {
		for (uint i = 0;i < rownum;++i)
			my_matrix[i].resize(new_colnum);
		colnum = new_colnum;
	}
	else if (new_colnum > colnum) {
		while (colnum < new_colnum)
			(*this).add_col(value);
	}
}

template<typename entry_element>
vector<entry_element> matrix<entry_element>::matrix_to_vector() const
{
	if ((*this).is_empty())
		return vector<entry_element>();
	vector<entry_element>result(rownum*colnum, entry_element());
	for (uint i = 0;i < rownum;++i)
		for (uint j = 0;j < colnum;++j)
			result[i*colnum + j] = my_matrix[i][j];
	return result;
}

//operator overload
template<typename entry_element>
template<typename T>
matrix<entry_element> matrix<entry_element>::operator+(const matrix<T>& m) const
{
	if (rownum != m.row_num() || colnum != m.col_num()) {
		cout << "Matrix addition error: two matrix has different size\n";
		return matrix<entry_element>();
	}
	matrix<entry_element> result = *this;
	for (uint i = 0;i < rownum;++i)
		for (uint j = 0;j < colnum;++j)
			result[i][j] += m[i][j];
	return result;
}

template<typename entry_element>
template<typename T>
matrix<entry_element> matrix<entry_element>::operator-(const matrix<T>& m) const
{
	if (rownum != m.row_num() || colnum != m.col_num()) {
		cout << "Matrix substitution error: two matrix has different size\n";
		return matrix<entry_element>();
	}
	matrix<entry_element> result = *this;
	for (uint i = 0;i < rownum;++i)
		for (uint j = 0;j < colnum;++j)
			result[i][j] -= m[i][j];
	return result;
}

template<typename entry_element>
template<typename T>
matrix<entry_element> matrix<entry_element>::operator*(const matrix<T>& m) const
{
	if (colnum != m.row_num()) {
		cout << "matrix multiplication error: invalid size\n";
		return matrix<entry_element>();
	}
	matrix<entry_element> result(rownum, m.col_num(), entry_element());
	for (uint i = 0;i < rownum;++i)
		for (uint j = 0;j < m.col_num();++j)
			for (uint k = 0;k < colnum;++k)
				result[i][j] += my_matrix[i][k] * m[k][j];
	return result;
}

template<typename entry_element>
template<typename T>
matrix<entry_element> matrix<entry_element>::operator*(const T & scalar) const
{
	matrix<entry_element> result = *this;
	for (uint i = 0;i < rownum;++i)
		for (uint j = 0;j < colnum;++j)
			result[i][j] *= scalar;
	return result;
}

template<typename entry_element>
template<typename T>
matrix<entry_element> matrix<entry_element>::operator/(const T & scalar) const
{
	matrix<entry_element> result = *this;
	for (uint i = 0;i < rownum;++i)
		for (uint j = 0;j < colnum;++j)
			result[i][j] /= scalar;
	return result;
}

template<typename entry_element>
matrix<double> matrix<entry_element>::operator%(const matrix<entry_element>& vect) const
{
	if (rownum != colnum) {
		cout << "operator% error: matrix is not square\n";
		return matrix<double>();
	}
	if (rownum != vect.row_num()) {
		cout<< "operator% error: matrix and vect has differect size\n";
		return matrix<double>();
	}
	matrix<double> aug_m(*this);//augmented matrix
	aug_m.resize(rownum, colnum + vect.col_num());
	for (uint i = 0;i < rownum;++i)
		for(uint j=0;j<vect.col_num();++j)
		aug_m[i][colnum+j] = double(vect[i][j]);
	aug_m.Gaussian_Diagonalization();
	return aug_m.submatrix(0, colnum, rownum, vect.col_num());
}

template<typename entry_element>
matrix<entry_element> matrix<entry_element>::operator-() const
{
	if ((*this).is_empty())
		return *this;
	matrix<entry_element> result(*this);
	for (uint i = 0;i < rownum;++i)
		for (uint j = 0;j < colnum;++j)
			result[i][j] = -result[i][j];
	return result;
}

template<typename entry_element>
matrix<entry_element> matrix<entry_element>::operator~() const
{
	if ((*this).is_empty())
		return *this;
	matrix<entry_element> result(colnum, rownum,entry_element());
	for (uint i = 0;i < rownum;++i)
		for (uint j = 0;j < colnum;++j)
			result[j][i] = my_matrix[i][j];
	return result;
}
//matrix operator

template<typename entry_element>
void matrix<entry_element>::add_row(const entry_element & value)
{
	if ((*this).is_empty())
		return;
	my_matrix.push_back(vector<entry_element>(colnum, value));
	++rownum;
}

template<typename entry_element>
void matrix<entry_element>::add_row(int location, const entry_element & value)
{
	if ((*this).is_empty())
		return;
	if (location<0 || location>rownum) {
		cout << "add_row error: invalid location\n";
		return;
	}
	if (location == rownum) { //add at the end
		(*this).add_row(value);
		return;
	}
	int loca_of_lastrow = rownum - 1; //index for the original last row
	my_matrix.push_back(my_matrix[loca_of_lastrow]); //copy the last row and add to the end
	int num_of_move = loca_of_lastrow - location;
	for (;num_of_move > 0;--num_of_move) { //move all the row between add location and the last row one unit down
		my_matrix[loca_of_lastrow] = my_matrix[loca_of_lastrow - 1];
		--loca_of_lastrow;
	}//loca_of_lastrow should point to the location of adding now
	for (int j = 0;j < colnum;++j) //change the value of that row
		my_matrix[loca_of_lastrow][j] = value;
	++rownum;
}

template<typename entry_element>
void matrix<entry_element>::add_row(const vector<entry_element>& v)
{
	if (v.empty())
		return;
	if ((*this).is_empty()) {
		my_matrix.push_back(v);
		rownum = 1;
		colnum = v.size();
		return;
	}
    if (v.size() > colnum)
		my_matrix.push_back(vector<entry_element>(&v[0], &v[colnum])); //truncation
	else {
		my_matrix.push_back(v);
		while (my_matrix[rownum].size() < colnum) //add zero 
			my_matrix[rownum].push_back(entry_element());
	}
	++rownum;
}

template<typename entry_element>
void matrix<entry_element>::add_row(int location, const vector<entry_element>& v)
{
	if (v.empty())
		return;
	if (location<0 || location>rownum) {
		cout << "add_row error: invalid location\n";
		return;
	}
	if (location == rownum) { //add at the end
		(*this).add_row(v);
		return;
	}
	int loca_of_lastrow = rownum - 1; //index for the original last row
	my_matrix.push_back(my_matrix[loca_of_lastrow]); //copy the last row and add to the end
	int num_of_move = loca_of_lastrow - location;
	for (;num_of_move > 0;--num_of_move) { //move all the row between add location and the last row one unit down
		my_matrix[loca_of_lastrow] = my_matrix[loca_of_lastrow - 1];
		--loca_of_lastrow;
	}//loca_of_lastrow should point to the location of adding now
	for (int j = 0;j < colnum;++j) //change the value of that row
		my_matrix[loca_of_lastrow][j] = j<v.size()?v[j]:entry_element();
	++rownum;
}

template<typename entry_element>
void matrix<entry_element>::remove_row(int location)
{
	if ((*this).is_empty())
		return;
	if (location<0 || location>=rownum) {
		cout << "remove_row error: invalid location\n";
		return;
	}
	if (rownum == 1) { //case of the matrix only has one row
		(*this).clear();
		return;
	}
	int loca_of_lastrow = rownum - 1;
	while (location != loca_of_lastrow) { //move the deleted row one unit down until it reach the bottom
		my_matrix[location] = my_matrix[location + 1];
		++location;
	}
	my_matrix.pop_back();
	--rownum;
}

template<typename entry_element>
void matrix<entry_element>::add_col(const entry_element & value)
{
	if ((*this).is_empty())
		return;
	for (uint i = 0;i < rownum;++i)
		my_matrix[i].push_back(value);
	++colnum;
}

template<typename entry_element>
void matrix<entry_element>::add_col(int location, const entry_element & value)
{
	if ((*this).is_empty())
		return;
	if (location<0 || location>colnum) {
		cout << "add_col error: invaild location\n";
		return;
	}

	(*this).add_col(value); //add the new col at the end first
	int cur_col = colnum-1; //current colnum's index
	if (cur_col == location)
		return;
	while (cur_col > location) {
		col_copy(cur_col, cur_col - 1);
		--cur_col;
	}
	for (uint i = 0;i < rownum;++i)
		my_matrix[i][cur_col] = value;
}

template<typename entry_element>
void matrix<entry_element>::add_col(const vector<entry_element>& v)
{
	if (v.empty())
		return;
	if ((*this).is_empty())
		(*this) = matrix<entry_element>(v);
	else {
		uint loca_of_lastcol = colnum;
		(*this).add_col();
		uint range = v.size() < rownum ? v.size() : rownum;
		for (uint i=0;i < range;++i)
			my_matrix[i][loca_of_lastcol] = v[i];
	}
}

template<typename entry_element>
void matrix<entry_element>::add_col(int location, const vector<entry_element>& v)
{
	if (v.empty())
		return;
	if (location<0 || location>colnum) {
		cout << "add_col error: invaild location\n";
		return;
	}
	if ((*this).is_empty())
		(*this).add_col(v);
	else {
		uint cur_col = colnum;
		(*this).add_col();
		while (cur_col > location) {
			col_copy(cur_col, cur_col - 1);
			--cur_col;
		}
		for (uint i = 0;i < rownum;++i)
			my_matrix[i][cur_col] = i<v.size()?v[i]:entry_element();
	}
}

template<typename entry_element>
void matrix<entry_element>::remove_col(int location)
{
	if ((*this).is_empty())
		return;
	if (location<0 || location>colnum - 1) {
		cout << "remove_col error: invaild location\n";
		return;
	}
	if (colnum == 1) {
		(*this).clear();
		return;
	}
	while (location < colnum - 1) {
		col_copy(location, location + 1);
		++location;
	}
	for (uint i = 0;i < rownum;++i)
		my_matrix[i].pop_back();
	--colnum;
}

template<typename entry_element>
void matrix<entry_element>::add_entry(int row_index, int col_index, const entry_element & value)
{
	if (row_index < 0 || col_index < 0) {
		cout << "add_entry error: invalid row_index and/or col_index\n";
		return;
	}
	if ((*this).is_empty())
		(*this).resize(row_index + 1, col_index + 1);
	while (rownum <= row_index)
		(*this).add_row();
	while (colnum <= col_index)
		(*this).add_col();
	my_matrix[row_index][col_index] = value;
}

template<typename entry_element>
matrix<entry_element> matrix<entry_element>::get_row(int row_index) const
{
	if (row_index < 0 || row_index >= rownum) {
		cout << "get_row error: invalid row index\n";
		return matrix<entry_element>();
	}
	matrix<entry_element> result(1,colnum,entry_element()); //a matrix with only 1 row
	for (uint j = 0;j < colnum;++j)
		result[0][j] = my_matrix[row_index][j];
	return result;
}

template<typename entry_element>
matrix<entry_element> matrix<entry_element>::get_col(int col_index) const
{
	if (col_index < 0 || col_index >= colnum) {
		cout << "get_col error: invalid column index\n";
		return matrix<entry_element>();
	}
	matrix<entry_element> result(rownum, 1, entry_element()); //a matrix with only 1 col
	for (uint i = 0;i < rownum;++i)
		result[i][0] = my_matrix[i][col_index];
	return result;
}

template<typename entry_element>
matrix<entry_element> matrix<entry_element>::submatrix(int row_index, int col_index, int num_of_row, int num_of_col)
{
	if (row_index < 0 || col_index < 0 || num_of_row < 0 || num_of_col < 0||row_index+num_of_row>rownum||col_index+num_of_col>colnum) {
		cout << "submatrix error: invalid input arguments\n";
		return matrix<entry_element>();
	}
	matrix<entry_element> result(num_of_row, num_of_col, entry_element());
	for (uint i = 0;i < num_of_row;++i)
		for (uint j = 0;j < num_of_col;++j)
			result[i][j] = my_matrix[row_index + i][col_index + j];
	return result;
}

//friend function
template<typename entry_element>
matrix<double> inv(const matrix<entry_element>& m)
{
	if (m.rownum != m.colnum) {
		cout << "matrix inv error: matrix is not square\n";
		return matrix<double>();
	}
	matrix<entry_element> aug_m(m); //augmented
	aug_m.resize(m.rownum, 2 * m.colnum, entry_element());
	for (uint i = 0;i < m.rownum;++i)
		aug_m[i][i + m.colnum] = entry_element(1);
	if (!aug_m.Gaussian_Diagonalization()) {
		cout<< "matrix inv error : matrix is singularly\n";
		return matrix<double>();
	}
	return matrix<double> (aug_m.submatrix(0, m.colnum, m.rownum, m.rownum));
}
template<typename T,typename entry_element>
matrix<entry_element> operator*(T scalar, const matrix<entry_element>& m)
{
	matrix<entry_element> result = m;
	for (uint i = 0;i < m.row_num();++i)
		for (uint j = 0;j < m.col_num();++j)
			result[i][j] *= scalar;
	return result;
}

//non-mumber function

template<typename entry_element>
matrix<entry_element> urand(int row_num, int col_num, const entry_element & begin, const entry_element & end)
{
	if (row_num < 0 || col_num < 0 || end < begin) {
		cout << "urand error: invalid input arguments\n";
		return matrix<entry_element>();
	}
	if (row_num == 0 || col_num == 0)
		return matrix<entry_element>();
	matrix<entry_element> result(row_num, col_num, entry_element());
	for(uint i=0;i<row_num;++i)
		for(uint j=0;j<col_num;++j)
			result[i][j]= (double(rand()) / RAND_MAX)*(end - begin) + begin;
	return result;
}

template<typename entry_element>
matrix<entry_element> reshape(const vector<entry_element>& vect, int row_num, int col_num)
{
	if (vect.empty())
		return matrix<entry_element>();
	if (vect.size() != row_num*col_num) {
		cout << "reshape error: incorrect element number\n";
		return matrix<entry_element>();
	}
	matrix<entry_element> result(row_num, col_num, entry_element());
	for (uint i = 0;i < row_num;++i)
		for (uint j = 0;j < col_num;++j)
			result[i][j] = vect[i*col_num + j];
	return result;
}
#endif
