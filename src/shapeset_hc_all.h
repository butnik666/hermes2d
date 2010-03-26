// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __HERMES2D_SHAPESET_HC_ALL_H
#define __HERMES2D_SHAPESET_HC_ALL_H

#ifdef COMPLEX

#include "shapeset.h"


/// H(curl) shapeset based on Legendre polynomials.
class HERMES2D_API HcurlShapesetLegendre : public Shapeset
{
  public: HcurlShapesetLegendre();
  virtual int get_id() const { return 10; }
};


// Experimental.
class HERMES2D_API HcurlShapesetEigen2 : public Shapeset
{
  public: HcurlShapesetEigen2();
  virtual int get_id() const { return 11; }
};


// Experimental.
class HERMES2D_API HcurlShapesetGradEigen : public Shapeset
{
  public: HcurlShapesetGradEigen();
  virtual int get_id() const { return 12; }
};


/// H(curl) shapeset with Legendre bubbles and gradients of H1 functions as edges
class HERMES2D_API HcurlShapesetGradLeg : public Shapeset
{
  public: HcurlShapesetGradLeg();
  virtual int get_id() const { return 13; }
};


/// This is the default Hcurl shapeset typedef.
typedef HcurlShapesetGradLeg HcurlShapeset;


#endif // COMPLEX

#endif
