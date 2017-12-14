#ifndef BG_LARRAY_H_
#define BG_LARRAY_H_
#include <sparsehash/dense_hash_map> // or sparse_hash_set, dense_hash_map, ...
#include <google/dense_hash_map>
#include "limits.h"
using google::dense_hash_map;

#define DENSE_HASHMAP

typedef int vertex;
template <class T>
class HashMap
{
	private:
		T initialValue_;
		T* array;
    	dense_hash_map<size_t, T> data_;

	public:
    	typedef typename dense_hash_map<size_t,T>::iterator iterator;
    	T & operator [] (size_t index) {
			return data_[index];
		}
        void setEmptyKey() {
        	data_.set_empty_key(-1);
        }
        void setDeletedKey() {
        	data_.set_deleted_key(-2);
        }
    	bool hasDefaultValue(size_t index) {
			iterator it = data_.find(index);
			if (it==data_.end())
				return true;
			return (it->second == initialValue_);
    	}
    	HashMap(T initialValue) :
    		initialValue_(initialValue) {
    		data_.set_empty_key(-1);
    	}
		void reset(T initialValue) {
			data_.clear();
			initialValue_ = initialValue;
		}
        void erase(size_t index) {
        		data_.erase(index);
        }
    	int size() {
    		return data_.size();
    	}
        iterator begin() {
            return data_.begin();
        }
        iterator end() {
            return data_.end();
        }
};

#endif
