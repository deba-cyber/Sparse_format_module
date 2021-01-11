
# include "sparse_module.hpp"
# include <iostream>
# include <string>
# include <vector>
# include <algorithm>



/****************************************************
 * Given rowptr/colptr , following routine will 
 * generate rowind/colind ..........
 * to be used to generate COO (rowind/colind)
 * from CSR/CSC rowptr/colptr 
 * **************************************************/

std::vector<int> rowind_or_colind_from_ptr(const std::vector<int>& inputvect)	{
	std::vector<int> tmpvect;
	for (int i=1; i<inputvect.size(); i++)	{
		int tmpcur = inputvect[i]-inputvect[i-1];
		tmpvect.insert(tmpvect.end(),tmpcur,i-1);
	}
	return tmpvect;
}

/****************************************************/
/****************************************************/

/* Implementing functions of class CSC_CSR_COO_BACK_N_FORTH */

////				Constructor		///////

CSC_CSR_COO_BACK_N_FORTH::CSC_CSR_COO_BACK_N_FORTH(const std::string& sparsetype)	{
	this-> sparsetype = sparsetype;
}

/*************************************************
 * member function to generate CSR structure from 
 * COO ...
 * input vectors .. 
 * first_ind_vect => row indices 
 * second_ind_vect => column indices 
 * for CSR generation column indices and data 
 * will be sorted in place .................
 * only new vector is row pointer vector to be
 * returned
 *************************************************/

std::vector<int> CSC_CSR_COO_BACK_N_FORTH::generate_CSR_from_COO(std::vector<int>& first_ind_vect,
		std::vector<int>& second_ind_vect,std::vector<double>& data_vect,const std::string& rtrn)	{
	/**************************************************
	 * when input type is COO => first_ind_vect &  *
	 * second_ind_vect are ROW & COLUMN indices 
	 * both are 0-based indexing ......................
	 * ************************************************/
	if (rtrn != "CSR"	||	sparsetype != "COO")	{
		throw "Check input or output sparse format types";
	}
	else	{
		sort_vectors_w_bubblesort(first_ind_vect,second_ind_vect,data_vect);
		std::vector<int> rowptr;
		int row_size = *max_element(first_ind_vect.begin(),first_ind_vect.end());
		rowptr.reserve(row_size+2);		// adding 2 as highest element has 0-based indexing
		rowptr.emplace_back(0);
		int cur_ind = first_ind_vect[0];
		int counter = 1;
		for (int i=1; i<first_ind_vect.size(); i++)	{
			if (first_ind_vect[i] != cur_ind && i!=first_ind_vect.size()-1)	{
				rowptr.emplace_back(counter);
				counter = 1;
				cur_ind = first_ind_vect[i];
			}
			else if (i == first_ind_vect.size()-1 && first_ind_vect[i] == first_ind_vect[i-1])	{
				counter += 1;
				rowptr.emplace_back(counter);
				break;
			}
			else if (i == first_ind_vect.size()-1 && first_ind_vect[i] != first_ind_vect[i-1])	{
				rowptr.emplace_back(counter);
				rowptr.emplace_back(1);
				break;
			}
			else	{
				counter += 1;
				cur_ind = first_ind_vect[i];
				continue;
			}
		}
		std::vector<int> rowptr_rtrn = partial_sum(rowptr);
		reorder_from_ptrvect(second_ind_vect,data_vect,rowptr_rtrn);
		return rowptr_rtrn;
	}
}

/***************************************************
 * member function to generate CSC structure from 
 * COO ...
 * input vectors ..
 * first_ind_vect => row indices 
 * second_ind_vect => column indices
 * for CSC generation row indices and data will be 
 * sorted in place ...............................
 * only new vector will be column pointer ........
 * **************************************************/

std::vector<int> CSC_CSR_COO_BACK_N_FORTH::generate_CSC_from_COO(std::vector<int>& first_ind_vect,
		std::vector<int>& second_ind_vect,std::vector<double>& data_vect,const std::string& rtrn)	{
	/**************************************************
	 * when input type is COO => first_ind_vect &  *
	 * second_ind_vect are ROW & COLUMN indices 
	 * both are 0-based indexing ......................
	 * ************************************************/
	if (rtrn != "CSC"	||	sparsetype != "COO")	{
		throw "Check input or output sparse format types";
	}
	else	{
		sort_vectors_w_bubblesort(second_ind_vect,first_ind_vect,data_vect);
		std::vector<int> colptr;
		int col_size = *max_element(second_ind_vect.begin(),second_ind_vect.end());
		colptr.reserve(col_size+2);		// adding 2 as highest element has 0-based indexing 
		colptr.emplace_back(0);
		int cur_ind = second_ind_vect[0];
		int counter = 1;
		for (int i=1; i<second_ind_vect.size(); i++)	{
			if (second_ind_vect[i] != cur_ind && i!=second_ind_vect.size()-1)	{
				colptr.emplace_back(counter);
				counter = 1;
				cur_ind = second_ind_vect[i];
			}
			else if (i == second_ind_vect.size()-1 && second_ind_vect[i] == second_ind_vect[i-1])	{
				counter += 1;
				colptr.emplace_back(counter);
				break;
			}
			else if (i == second_ind_vect.size()-1 && second_ind_vect[i] != second_ind_vect[i-1])	{
				colptr.emplace_back(counter);
				colptr.emplace_back(1);
				break;
			}
			else	{
				counter += 1;
				cur_ind = second_ind_vect[i];
				continue;
			}
		}
		std::vector<int> colptr_rtrn = partial_sum(colptr);
		reorder_from_ptrvect(first_ind_vect,data_vect,colptr_rtrn);
		return colptr_rtrn;
	}
}


/**********************************************************
 * member function to generate row_ind vector (sorted) 
 * from passed rowptr vector for generating CSR format
 * ********************************************************/


std::vector<int> CSC_CSR_COO_BACK_N_FORTH::generate_COO_from_CSR(std::vector<int>& rowptr_vect,const std::string& rtrn)	{
	if (rtrn != "COO" || sparsetype != "CSR")	{
		throw	"Check input or output sparse format types";
	}
	else	{
		std::vector<int> row_ind_vect = rowind_or_colind_from_ptr(rowptr_vect);
		return row_ind_vect;
	}
}


/**********************************************************
 * member function to generate col_ind vector (sorted)
 * from passed colptr vector for generating CSC format
 * ********************************************************/

std::vector<int> CSC_CSR_COO_BACK_N_FORTH::generate_COO_from_CSC(std::vector<int>& colptr_vect,const std::string& rtrn)	{
	/**********************************************************
	 * ********************************************************/
	if (rtrn != "COO"	|| sparsetype != "CSC")	{
		throw	"Check input or output sparse format types";
	}
	else	{
		std::vector<int> col_ind_vect = rowind_or_colind_from_ptr(colptr_vect);
		return col_ind_vect;
	}
}




