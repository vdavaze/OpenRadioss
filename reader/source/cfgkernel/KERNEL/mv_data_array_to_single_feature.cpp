/*Copyright>    OpenRadioss
//Copyright>    Copyright (C) 1986-2025 Altair Engineering Inc.
//Copyright>
//Copyright>    This program is free software: you can redistribute it and/or modify
//Copyright>    it under the terms of the GNU Affero General Public License as published by
//Copyright>    the Free Software Foundation, either version 3 of the License, or
//Copyright>    (at your option) any later version.
//Copyright>
//Copyright>    This program is distributed in the hope that it will be useful,
//Copyright>    but WITHOUT ANY WARRANTY; without even the implied warranty of
//Copyright>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//Copyright>    GNU Affero General Public License for more details.
//Copyright>
//Copyright>    You should have received a copy of the GNU Affero General Public License
//Copyright>    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//Copyright>
//Copyright>
//Copyright>    Commercial Alternative: Altair Radioss Software
//Copyright>
//Copyright>    As an alternative to this open-source version, Altair also offers Altair Radioss
//Copyright>    software under a commercial license.  Contact Altair to discuss further if the
//Copyright>    commercial version may interest you: https://www.altair.com/radioss/.*/

#include <UTILS/win32_utils.h>  

#include "mv_descriptor.h"
#include "mv_data_array_to_single_feature.h"


/* --------- Constructors & destructor --------- */

MvDataArrayToSingleFeature_t::MvDataArrayToSingleFeature_t(const string &name,int nb_features,int index) :
  MvDataArrayFeature_t(DFT_ARRAY_TO_SINGLE,name,nb_features,0),
  myIndex(index)
{}

MvDataArrayToSingleFeature_t::~MvDataArrayToSingleFeature_t() {
}


/* --------- Output in an output stream --------- */

ostream &MvDataArrayToSingleFeature_t::display(ostream &os,const MvDescriptor_t &descr,int level) const {
  int i;
  for(i=0;i<level;i++) os << "  ";
  display_props(os); 
  //
  os << "ARRAY_TO_SINGLE(TITLE=\"" << getTitle() << "\""
     << ",INDEX="                  << getIndex() 
     << ") {"                      << endl;
  //
  bool a_is_optional=false;
  for(i=0;i<getNumber();++i) {
    const MvDataFeature_t *a_df_p=getDataFeature(i);
    if(i==0 || a_df_p->isOptional()!=a_is_optional) {
      a_is_optional=a_df_p->isOptional();
      for(int j=0;j<level;j++) os << "  ";
      if(a_is_optional) os << "optional:"; else os << "mandatory:";
      os << endl;
    }
    a_df_p->display(os,descr,level+1) << endl;
  }
  //
  for(i=0;i<level;i++) os << "  ";
  os << "}";
  return os;
}




