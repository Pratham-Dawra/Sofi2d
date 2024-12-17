
/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI2D.
 * 
 * SOFI2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * SOFI2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SOFI2D. See file COPYING and/or 
  * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 * Initialization of the wave field with zero values (zero wavefield)
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void zero_acoustic(int j, int i, MemWavefield * mpw)
{
    mpw->pvx[j][i] = 0.0f;
    mpw->pvy[j][i] = 0.0f;
    mpw->psxx[j][i] = 0.0f;
    mpw->psyy[j][i] = 0.0f;
}
