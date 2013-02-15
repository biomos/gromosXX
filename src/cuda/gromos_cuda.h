/*
 * File:   gromos_cuda.h
 * Author: Axel && Mel
 *
 * Created on October 20, 2011, 11:34 AM
 *
 * Global CUDA include file
 */

#ifndef GROMOS_CUDA_H
#define	GROMOS_CUDA_H
#include <memory>
#include <new>
#include <iostream>
#include <map>
//#include <cstdint>

// define namespace
namespace gcuda {

    // Heap allocation class for a 2D array on the host
    // This allacates (width*height) of linear memory
    // in order to avoid the costly new as much as possible, a one dimentional array of size
    // width is created to hold the accession pointers. Then another one dimentional
    // array of size (width*height) is created at position [0] of the first array and
    // the first array positions 1 to width are filled with pointers to the first array
    // with [0] + i * height

    template <class T>
    class TwoDArray {
        public:
            // Accession pointer
            T** ptr;
            // Memory pitch
            size_t pitch;
            // Memory allocation error
            bool error;
            // Size array
            // This stores how many dim_j elements acutally hold values
            unsigned int* elements;

            // the dimensions of the array
            unsigned int dim_i;
            unsigned int dim_j;

            /*
             * Constructor
             */
            TwoDArray(unsigned int width, unsigned int height){
            error = false;
            dim_i = width;
            dim_j = height;
            //const std::nothrow_t nothrow;
            ptr = new(std::nothrow) T* [width];
                if (ptr == NULL) {
                        error = true;
                }
                ptr[0] = new(std::nothrow) T [width*height];
                if (ptr[0] == NULL) {
                        error = true;
                } else {
                        for (unsigned int i = 1; i < width; i++) {
                                ptr[i] = ptr[0] + i * height;

                        }
                        pitch = height * sizeof(T);
                }
                elements = new(std::nothrow) unsigned int [width];
	    	if (elements == NULL) {
                	error = true;
                }
            }

            /*
             * Deconstructor
             */
            ~TwoDArray(){
            delete[] ptr[0];
            delete[] ptr;
            delete[] elements;
            }
    };

    // iterator class for the TwoDArray class
    template <class T>
    class TwoDArray_iterator : public std::iterator<std::forward_iterator_tag, T>
{
public:
  TwoDArray_iterator() {}
  TwoDArray_iterator(unsigned int* i) {
      it = i;
  }
  TwoDArray_iterator(unsigned int** i) {
      it = *i;
  }
  TwoDArray_iterator(const TwoDArray_iterator& tda_it) {
      it = tda_it.it;
  }
  TwoDArray_iterator& operator++() {
      ++it;
      return *this;
  }
  TwoDArray_iterator operator++(int) {
        TwoDArray_iterator tmp(*this);
        operator++();
        return tmp;
  }
  TwoDArray_iterator& operator+=(int i) {
      TwoDArray_iterator tmp(*this);
      for (unsigned int j = 0; j < i; j++){
        ++tmp;
      }
      return tmp;
  }
  bool operator==(const TwoDArray_iterator& tda_it) {
      return it==tda_it.it;
  }
  bool operator!=(const TwoDArray_iterator& tda_it) {
      return it!=tda_it.it;
  }
  unsigned int& operator*() {
      return *it;
  }
    private:
        T* it;
};

// memory map class, this holds the pointers to the device memory - sort of a map
// of the gpu memory, it uses std::map to point to an array of pointers, the array is
// needed to be abled to be copied to the device

class memmap {

public:
    std::map<std::string, unsigned int> map;
    //uintptr_t* memory;
    size_t* memory;
    unsigned int size;
    unsigned int counter;
    size_t* pitches;

    memmap(unsigned int elements) {
        //should be malloc
        //change with new compiler
        //memory = new(std::nothrow) uintptr_t [elements];
        memory = new(std::nothrow) size_t [elements];
        pitches = new(std::nothrow) size_t [elements];
        size = elements;
        counter = 0;
    }
    
    ~memmap() {
        delete[] memory;
        delete[] pitches;
    }

    void add(std::string description, void* pointer, size_t pitch) {
        map.insert(std::pair<std::string,unsigned int>(description, counter));
        //memory[counter] = (uintptr_t)pointer;
        memory[counter] = (size_t)pointer;
        pitches[counter] = pitch;
        counter++;
    }
};

}//namespace


#endif	/* GROMOS_CUDA_H */