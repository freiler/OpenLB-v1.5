/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Mathias J. Krause, Cyril Masquelier,
 *  Benjamin Förster, Albert Mink
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

#ifndef SMOOTH_INDICATOR_CALC_F_2D_HH
#define SMOOTH_INDICATOR_CALC_F_2D_HH

#include "smoothIndicatorCalcF2D.h"

namespace olb {


//////////////////////////////// IndicSmoothCalc2D ////////////////////////////////
template <typename T, typename S>
SmoothIndicCalc2D<T,S>::SmoothIndicCalc2D(SmoothIndicatorF2D<T,S>& f, SmoothIndicatorF2D<T,S>& g)
  : _f(f), _g(g)
{
  for ( int i=0; i<2; i++) {
    this->_myMin[i] = util::min(f.getMin()[i], g.getMin()[i]);
    this->_myMax[i] = util::max(f.getMax()[i], g.getMax()[i]);
  }
  this->_circumRadius = util::max( 0.5*(this->getMax()[0]-this->getMin()[0]), 0.5*(this->getMax()[1]-this->getMin()[1]));
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}


template <typename T, typename S>
SmoothIndicPlus2D<T,S>::SmoothIndicPlus2D(SmoothIndicatorF2D<T,S>& f, SmoothIndicatorF2D<T,S>& g)
  : SmoothIndicCalc2D<T,S>(f,g)
{}

// returns 1 if( f==1 || g==1 ) UNION
template <typename T, typename S>
bool SmoothIndicPlus2D<T, S>::operator()(T output[], const S input[])
{
  this->_f(output, input);
  T tmp;
  this->_g(&tmp, input);
  output[0] = util::max(output[0], tmp);
  return true;
}


template <typename T, typename S>
SmoothIndicatorF2D<T,S>& SmoothIndicatorF2D<T,S>::operator+(SmoothIndicatorF2D<T,S>& rhs)
{
  auto tmp = std::make_shared< SmoothIndicPlus2D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}


} // end namespace olb

#endif
