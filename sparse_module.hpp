# if !defined (SPARSE_MODULE)
# define SPARSE_MODULE

# include <iosfwd>
# include <string>
# include <vector>
# include <algorithm>



/*********************************************
 * template function to sort row/col indices
 * using bubblesort .. 
 * other two vectors col/row & data will be 
 * sorted in the same order in place
 * this will generate rowind/colind .....
 * still to be sorted for each row/column
 * written in a way such that sorting can 
 * be done when first two vectors are of same 
 * data type and third one has a different data 
 * type	.......................................
 * ********************************************/

template<typename T1,typename T2>
inline void sort_vectors_w_bubblesort(std::vector<T1>& firstvect, 
		std::vector<T1>& followvect1, std::vector<T2>& followvect2)	{
	int vectsize = firstvect.size();
	std::string swapped;
	for (int i=0; i<vectsize-1; i++)	{
		swapped = "False";
		for (int j=0; j<vectsize-i-1; j++)	{
			if (firstvect[j] > firstvect[j+1])	{
				std::swap(firstvect[j],firstvect[j+1]);
				std::swap(followvect1[j],followvect1[j+1]);
				std::swap(followvect2[j],followvect2[j+1]);
				swapped = "True";
			}
		}
		if (swapped == "False")	
			break;
	}
}

/**************************************************
 * template function to do partial sum in a vector
 * generic ... 
 * here to be used to generate rowptr/colptr from 
 * vector containing no of non-zeros in each 
 * row/col .......................................
 * ***********************************************/

template<typename T>
inline std::vector<T> partial_sum(std::vector<T>& inputvect)	{
	std::vector<T> rtrnvect;
	rtrnvect.reserve(inputvect.size());
	rtrnvect.emplace_back(inputvect[0]);
	T cur_sum = 0.;
	for (int i=1; i<inputvect.size(); i++)	{
		cur_sum += inputvect[i];
		rtrnvect.emplace_back(cur_sum);
	}
	return rtrnvect;
}

/****************************************************
 * Given rowptr/colptr , following routine will 
 * generate rowind/colind ..........
 * to be used to generate COO (rowind/colind)
 * from CSR/CSC rowptr/colptr 
 * **************************************************/

std::vector<int> rowind_or_colind_from_ptr(const std::vector<int>& inputvect);

/*******************************************************
 * template function to generate final (colind & data)
 * i.e. sorted in each row for CSR 
 * and (rowind & data)	i.e. sorted in each column for 
 * CSC with passed inputs of colind/rowind with data 
 * (not sorted in each col/row) generated from bubblesort
 * ****************************************************/

template<typename T>
inline void reorder_from_ptrvect(std::vector<int>& row_or_col_ind_wo_sort,
		std::vector<T>& datavect_wo_sort,const std::vector<int>& row_or_col_ptr)	{
	std::vector<std::pair<int,T>> tmpvec;
	int size = datavect_wo_sort.size();
	for (int i=0; i<size; i++)	{
		tmpvec.push_back(std::make_pair(row_or_col_ind_wo_sort[i],datavect_wo_sort[i]));
	}
	for (int i=0; i<row_or_col_ptr.size()-1; i++)	{
		int start_id = row_or_col_ptr[i];
		int end_id = row_or_col_ptr[i+1];
		std::sort(tmpvec.begin()+start_id, tmpvec.begin()+end_id);
	}
	for (int i=0; i<size; i++)	{
		row_or_col_ind_wo_sort[i] = tmpvec[i].first;
		datavect_wo_sort[i] = tmpvec[i].second;
	}
}




/***************************************************
 * class for interconversion of sparse 
 * matrix structures CSR <-> CSC <-> COO 
 * string sparsetype gives input sparse type 
 * for all member functions for interconversion
 * three vectors are passed by reference ..
 * *************************************************/


class CSC_CSR_COO_BACK_N_FORTH	{
	private:
		std::string sparsetype;
	public:
		CSC_CSR_COO_BACK_N_FORTH(const std::string& sparsetype);	
		// ======== Member Functions ======== //
		std::vector<int> generate_CSR_from_COO(std::vector<int>& first_ind_vect,
				std::vector<int>& second_ind_vect,std::vector<double>& data_vect, const std::string& rtrn);
		std::vector<int> generate_CSC_from_COO(std::vector<int>& first_ind_vect,
				std::vector<int>& second_ind_vect,std::vector<double>& data_vect, const std::string& rtrn);
		std::vector<int> generate_COO_from_CSR(std::vector<int>& rowptr_vect,
				const std::string& rtrn);
		std::vector<int> generate_COO_from_CSC(std::vector<int>& colptr_vect,
				const std::string& rtrn);
};

#endif


