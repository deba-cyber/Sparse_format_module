# if ! defined (FILEOPS)
# define FILEOPS

# include <iostream>
# include <fstream>
# include <string>
# include <vector>
# include <algorithm>



/*****************************************************************
 * Function for accumulating data from multiple files with many
 * lines and single/multiple columns .. returning a single vector	
 * all values from all files in a single vector	.................
 * This function is specific to this calculation ................
 * working with the files containing transformation matrices 
 * SO => ANH & ANH_2_PO and prim_dvr_points files ...............
 * ***************************************************************/

std::vector<double> vect_accumulate_all_col(int* mode_id_arr, const std::string s1, int NDIM);

/////////////////////////////////////////////////////

/****************************************************
 * Function template to free up memory allocated 
 * for vector of any data type	********************
 * **************************************************/


template<typename T>
inline auto empty_swap(std::vector<T> & vec)	{
	std::vector<T>().swap(vec);
	return vec;
}

/////////////////////////////////////////////////////

/****************************************************
 * Function template to return specific queried row *
 * or column elements from passed one dimensional 
 * vector (which is obtained from two dimensional 
 * file data ........................***************
 * **************************************************/

template<typename T>
inline auto get_specific_row_or_col_elements(std::vector<T> & vec_onedim,
		int N_row_or_col, int N_tot_rows, int N_tot_cols, 
		const std::string s1)	
{
	/** Indexing for N_row_or_col is 0-based i.e. if 
	 *	N_row_or_col is 100 , actually 101-th row is desired **/
	std::vector<T> rtrn_sp_row_or_col_elems;
	std::vector<int> col_iterator_vect;		// index vector over which loop will run for column
	int start_cnt;
	int end_cnt;
	int cur_row = 0;
	int cur_id_4_col = N_row_or_col;		// index for iterating over for extracting column 
	if (s1 == "ROW")	{
		start_cnt = N_row_or_col*N_tot_cols;
		end_cnt = (N_row_or_col+1)*N_tot_cols;
		for (int i=start_cnt; i<end_cnt; i++)	{
			rtrn_sp_row_or_col_elems.emplace_back(vec_onedim[i]);
		}
		return rtrn_sp_row_or_col_elems;
	}
	else if (s1 == "COLUMN")	{
		while (cur_row != N_tot_rows)	{
			rtrn_sp_row_or_col_elems.emplace_back(vec_onedim[cur_id_4_col]);
			cur_id_4_col += N_tot_cols;
			cur_row +=1;
		}
		return rtrn_sp_row_or_col_elems;
	}
}

/*******************************************************
 * Function template to return matrix-vector product ***
 * where the matrix is passed in argument as 1-dim 
 * vector	by reference .............................
 * No of rows and columns are also passed by values 
 * *****************************************************/

template<typename T>
inline std::vector<T> get_mat_vec_pdt(std::vector<T> & mat_as_vect, std::vector<T> & onedim_vect,
		int N_row, int N_col)
{
	std::vector<T> rtrn_matvec_pdt_vect;
	std::vector<T> tmp_process_vect;
	int onedim_vect_size = onedim_vect.size();
	T sum;
	if (N_col != onedim_vect_size)	{
		throw "Matrix-Vector multiplication not possible!";
	}
	else
	{
		for (int i=0; i<N_row; i++)	{
			tmp_process_vect = get_specific_row_or_col_elements(mat_as_vect,i,N_row,N_col,"ROW");
			sum = 0.;
			for (int j=0; j<tmp_process_vect.size(); j++)	{
				sum += tmp_process_vect[j]*onedim_vect[j];
			}
			rtrn_matvec_pdt_vect.emplace_back(sum);
		}
	}
	return rtrn_matvec_pdt_vect;
}

/************************************************
 * Function for extracting a specific line from 
 * a file of arbitrary length	.. .............
 * this will be used to get the roots and weights
 * for GAUSS-HERMITE quadrature saved in file 
 * calculated before using numpy.polynomial routine
 * specific to data type .. here double
 * **********************************************/

std::vector<double> get_specific_line(int lineno, const std::string s1);	

/**************************************************
 * Function for extracting specific lines
 * from a file	.................................
 * string s1 is passed by value filename wo ext.
 * vector lineno_vect contains numbers in increasing
 * order corrsponsing to lines to be extracted	
 * lineno_vect has 1-based indexing	...............
 * returning as one dimensional vector ..........
 * specific to data type .. here double
 * ************************************************/

std::vector<double> get_all_sp_lines(std::vector<int> & lineno_vect, const std::string s1);	

/***************************************************************************
 * Function for returning stride array as a vector for any given vector as 
 * input argument which gives number of elements for each dimension ....
 * *************************************************************************/

std::vector<int> GET_STRIDE_ARR_4_ANY(const std::vector<int> & size_vect);	

/********************************************************************
 * template function for appending elements to an existing binary 
 * file ... 
 * data is passed by reference in an one dimensional vector
 * alongwith filename & mode which has to be "append" for this 
 * to work
 * ******************************************************************/

template<typename T>
inline void save_binary(std::vector<T> & vect,
		const std::string  filename, const std::string mode)	{
	std::string ext = ".bin";
	if (mode == "append")	{
		std::ofstream outfile(filename+ext, std::ios::binary|std::ios::app);
		if (outfile.is_open())	{
			outfile.write(reinterpret_cast<char*>(&vect[0]),vect.size()*sizeof(T));
		}
		else	{
			throw "Error opening file for writing!\n";
			exit(1);
		}	
	}
}

/***********************************************************************
 * template function for extracting specific elements from a binary 
 * file ....
 * start_elem_no gives the element from which extraction starts
 * (1-based indexing)
 * vectsize gives number of elements to be extracted including the 
 * starting element ....................................................
 * *********************************************************************/


template<typename T>
inline std::vector<T> rtrn_vec_from_bin(std::vector<T> & dummyvec,
		const std::string filename,int start_elem_no,
		int vectsize)	{
	/*****************************************
	 * passing a dummy vector of same data type 
	 * of zero size for keeping template structure
	 * ***************************************/
	std::ifstream infile;
	std::string ext = ".bin";
	infile.open(filename+ext, std::ios::binary);
	std::vector<T> vect(vectsize);
	infile.seekg((start_elem_no-1)*sizeof(T),std::ios::beg);
	if (infile.is_open())	{
		infile.read(reinterpret_cast<char*>(&vect[0]),vectsize*sizeof(T));
	}
	else	{
		throw "Error opening file for writing!\n";
		exit(1);
	}
	return vect;
}

/***************************************************
 * template function to get number of elements in 
 * a binary file ..................................
 * *************************************************/


template<typename T>
inline int get_num_elems(std::vector<T> & dummyvec,
		const std::string filename)	{
	/*****************************************
	 * passing a dummy vector of same data type 
	 * of size zero for keeping template structure
	 * ***************************************/
	std::ifstream infile;
	std::string ext = ".bin";
	infile.open(filename+ext,std::ios::binary);
	infile.seekg(0,std::ios::end);
	int filesize = infile.tellg();
	int num_elems = filesize/sizeof(T);
	return num_elems;
}


# endif
