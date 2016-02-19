
#include <vector>
#include <unordered_map>

template <typename ValueType>
class dense_array {
public:   
    typedef size_t SizeType;
    ValueType *vals;
    SizeType length;
    ValueType sentinal;
    
    dense_array(ValueType *array, SizeType size, ValueType sent) 
    : vals(array), length(size), sentinal(sent)
    {}
    
    template <typename IndexType>
    bool contains(const IndexType& index) const {
        if (index >= 0 && index < length) {
            return (vals[index]!=sentinal);
        } else { 
            return false;
        }
    }
    
    template <typename IndexType>
    SizeType count(const IndexType& index) const {
        return (SizeType)contains(index);
    }
    
    template <typename IndexType>
    ValueType get(const IndexType& index) const {
        return vals[index];
    }
    
    template <typename IndexType>
    void put(const IndexType& index, const ValueType& val) {
        vals[index] = val;
    }
};

template <typename ValueType> 
dense_array<ValueType> make_dense_array(ValueType* array, size_t length, 
                            ValueType sentinal) 
{
    return dense_array<ValueType>(array, length, sentinal);
}

/** This maps a vector with additional functions where a sentinal indicates 
 * absence from the array.
 */
template <typename DenseContainerType>
class dense_container {
public:   
    typedef size_t SizeType;
    typedef typename DenseContainerType::value_type ValueType;
    DenseContainerType container;
    ValueType sentinal;
    
    dense_container(DenseContainerType& array, ValueType sent) 
    : container(array),  sentinal(sent)
    {}
    
    template <typename IndexType>
    bool contains(const IndexType& index) const {
        if (index >= 0 && index < container.size()) {
            return (container[index]!=sentinal);
        } else { 
            return false;
        }
    }
    
    template <typename IndexType>
    SizeType count(const IndexType& index) const {
        return (SizeType)contains(index);
    }
    
    template <typename IndexType>
    ValueType get(const IndexType& index) const {
        return container[index];
    }
    
    template <typename IndexType>
    void put(const IndexType& index, const ValueType& val) {
        container[index] = val;
    }
};

template <typename ValueType> 
dense_container<ValueType> make_dense_container(ValueType* array, size_t length, 
                            ValueType sentinal) 
{
    return dense_array<ValueType>(array, length, sentinal);
}

template <typename SparseType>
struct const_sparse_array {
    const SparseType& map;
    typedef typename SparseType::data_type ValueType;
    ValueType sentinal;
    
    const_sparse_array(const SparseType& m, ValueType sent) 
    : map(m), sentinal(sent)
    {}
    
    template <typename IndexType>
    size_t count(const IndexType& index) const {
        return map.count(index);
    }
    
    template <typename IndexType>
    bool contains(const IndexType& index) const {
        if (count(index) > 0) {
            return true;
        } else { 
            return false;
        }
    }
    
    template <typename IndexType>
    ValueType get(const IndexType& index) const {
        typename SparseType::const_iterator got = map.find(index);
        if (got == map.end()) {
            return sentinal;
        } else {
            return got->second;
        }
    }
};

template <typename DenseType> 
dense_container<DenseType> make_dense_container(const DenseType& array, 
            typename DenseType::value_type sentinal) 
{
    return dense_container<DenseType>(array, sentinal);
}

template <typename SparseType> 
const_sparse_array<SparseType> make_const_sparse_array(const SparseType& array, 
            typename SparseType::data_type sentinal) 
{
    return const_sparse_array<SparseType>(array, sentinal);
}

template <typename SparseType>
struct sparse_array {
    SparseType& map;
    typedef typename SparseType::data_type ValueType;
    ValueType sentinal;
    
    sparse_array(SparseType& m, ValueType sent) 
    : map(m), sentinal(sent)
    {}
    
    template <typename IndexType>
    size_t count(const IndexType& index) const {
        return map.count(index);
    }
    
    template <typename IndexType>
    bool contains(const IndexType& index) const {
        if (count(index) > 0) {
            return true;
        } else { 
            return false;
        }
    }
    
    template <typename IndexType>
    ValueType get(const IndexType& index) const {
        typename SparseType::const_iterator got = map.find(index);
        if (got == map.end()) {
            return sentinal;
        } else {
            return got->second;
        }
    }
    
    template <typename IndexType>
    ValueType put(const IndexType& index, const ValueType& value) {
        map[index] = value;
    }
};

template <typename SparseType> 
sparse_array<SparseType> make_sparse_array(const SparseType& array, 
            typename SparseType::data_type sentinal) 
{
    return sparse_array<SparseType>(array, sentinal);
}

