/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013-2017 Albert Mink, Lukas Baron, Mathias J. Krause,
 *                          Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/
#ifndef BLOCK_CALC_F_2D_HH
#define BLOCK_CALC_F_2D_HH

#include "blockCalcF2D.h"

namespace olb {


template <typename T, template<typename> class F>
BlockCalcF2D<T,F>::BlockCalcF2D (BlockF2D<T>& f, BlockF2D<T>& g)
  : BlockF2D<T>(
      g.getBlockStructure(),
      f.getTargetDim() > g.getTargetDim() ? f.getTargetDim() : g.getTargetDim()),
    _f(f), _g(g),
    _glob{}, _fIsBlock(true), _gIsBlock(true)
{
  this->getName() = "(" + f.getName() + F<T>::symbol + g.getName() + ")";
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T, template<typename> class F>
BlockCalcF2D<T,F>::BlockCalcF2D (BlockF2D<T>& f, GenericF<T,int>& g, int glob)
  : BlockF2D<T>(
      f.getBlockStructure(),
      f.getTargetDim() > g.getTargetDim() ? f.getTargetDim() : g.getTargetDim()),
    _f(f), _g(g),
    _glob(glob), _fIsBlock(true), _gIsBlock(false)
{
  this->getName() = "(" + f.getName() + F<T>::symbol + g.getName() + ")";
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T, template<typename> class F>
BlockCalcF2D<T,F>::BlockCalcF2D (GenericF<T,int>& f, int glob, BlockF2D<T>& g)
  : BlockF2D<T>(
      g.getBlockStructure(),
      f.getTargetDim() > g.getTargetDim() ? f.getTargetDim() : g.getTargetDim()),
    _f(f), _g(g),
    _glob(glob), _fIsBlock(false), _gIsBlock(true)
{
  this->getName() = "(" + f.getName() + F<T>::symbol + g.getName() + ")";
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T, template<typename> class F>
bool BlockCalcF2D<T,F>::operator()(T output[], const int input[])
{
  T* outputF = output;
  T outputG[this->getTargetDim()];

  if ( this->_fIsBlock && this->_gIsBlock ) {
    this->_f(outputF, input);
    this->_g(outputG, input);
  }
  else {
    const int superInput[3] = { this->_glob, input[0], input[1] };

    if ( this->_fIsBlock ) {
      this->_f(outputF, input);
      this->_g(outputG, superInput);
    }
    else {
      this->_f(outputF, superInput);
      this->_g(outputG, input);
    }
  }

  if ( _f.getTargetDim() == 1 || _g.getTargetDim() == 1 ) {
    // scalar operation
    if ( _f.getTargetDim() == 1 ) {
      // apply the scalar f to possibly multidimensional g
      for (int i = 1; i < this->getTargetDim(); i++) {
        outputF[i] = outputF[0];
      }
    }
    else if ( _g.getTargetDim() == 1 ) {
      // apply scalar g to possibly multidimensional f
      for (int i = 1; i < this->getTargetDim(); i++) {
        outputG[i] = outputG[0];
      }
    }
  }

  for (int i = 0; i < this->getTargetDim(); i++) {
    output[i] = F<T>()(outputF[i], outputG[i]);
  }

  return true;
}


/////////////////////////////////operator()/// ////////////////////////////////
template <typename T>
BlockF2D<T>& BlockF2D<T>::operator+(BlockF2D<T>& rhs)
{
  auto tmp = std::make_shared< BlockCalcPlus2D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T>
BlockF2D<T>& BlockF2D<T>::operator-(BlockF2D<T>& rhs)
{
  auto tmp = std::make_shared< BlockCalcMinus2D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T>
BlockF2D<T>& BlockF2D<T>::operator*(BlockF2D<T>& rhs)
{
  auto tmp = std::make_shared< BlockCalcMultiplication2D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T>
BlockF2D<T>& BlockF2D<T>::operator/(BlockF2D<T>& rhs)
{
  auto tmp = std::make_shared< BlockCalcDivision2D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}


} // end namespace olb

#endif
