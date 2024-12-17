
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

void zero_wavefield(int iter, MemWavefield *mpw, MemInv * minv, GlobVar *gv, GlobVarInv *vinv)
{
    size_t nbyte_2 = (gv->NY + 2 * gv->ND) * (gv->NX + 2 * gv->ND) * sizeof(float);
    size_t nbyte_3 = nbyte_2 * gv->L;
    size_t nbyte_pml_x = gv->NY * 2 * gv->FW * sizeof(float);
    size_t nbyte_pml_y = gv->NX * 2 * gv->FW * sizeof(float);
    size_t nbyte_grad = gv->NY * gv->NX * sizeof(float);

    bzero(&(mpw->pvx[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
    bzero(&(mpw->pvy[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
    bzero(&(mpw->psxx[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
    bzero(&(mpw->psyy[-gv->ND + 1][-gv->ND + 1]), nbyte_2);

    if (gv->WEQ >= EL_ISO && gv->WEQ <= VEL_TTI) {  /*elastic cases */
        bzero(&(mpw->psxy[-gv->ND + 1][-gv->ND + 1]), nbyte_2);

        if (gv->MODE == FWI) {
            bzero(&(minv->ux[-gv->ND + 1][-gv->ND + 1]), nbyte_2); //added
            bzero(&(minv->uy[-gv->ND + 1][-gv->ND + 1]), nbyte_2); //added
            bzero(&(minv->uxy[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
            bzero(&(minv->uyx[-gv->ND + 1][-gv->ND + 1]), nbyte_2); //added
            bzero(&(minv->pvxp1[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
            bzero(&(minv->pvyp1[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
            bzero(&(minv->pvxm1[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
            bzero(&(minv->pvym1[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        
            /* initialize gradient matrices for each shot with zeros */
            bzero(&(minv->waveconv_shot[1][1]), nbyte_grad);
            bzero(&(minv->waveconv_rho_shot[1][1]), nbyte_grad);
            bzero(&(minv->waveconv_u_shot[1][1]), nbyte_grad);

            if ((vinv->EPRECOND == 1) || ((vinv->EPRECOND == 3) && (vinv->EPRECOND_ITER == iter || (vinv->EPRECOND_ITER == 0)))) {
                bzero(&(minv->Ws[1][1]), nbyte_grad);
                bzero(&(minv->Wr[1][1]), nbyte_grad);
                bzero(&(minv->We[1][1]), nbyte_grad);
            }
        }

        /* viscoelastic */
        if (gv->L) {
            bzero(&(mpw->pr[-gv->ND + 1][-gv->ND + 1][1]), nbyte_3);
            bzero(&(mpw->pp[-gv->ND + 1][-gv->ND + 1][1]), nbyte_3);
            bzero(&(mpw->pq[-gv->ND + 1][-gv->ND + 1][1]), nbyte_3);
        }

        /* PML Boundary */
        if (gv->ABS_TYPE == 1) {
            bzero(&(mpw->psi_sxx_x[1][1]), nbyte_pml_x);
            bzero(&(mpw->psi_sxy_x[1][1]), nbyte_pml_x);
            bzero(&(mpw->psi_vxx[1][1]), nbyte_pml_x);
            bzero(&(mpw->psi_vyx[1][1]), nbyte_pml_x);
            bzero(&(mpw->psi_vxxs[1][1]), nbyte_pml_x);

            bzero(&(mpw->psi_syy_y[1][1]), nbyte_pml_y);
            bzero(&(mpw->psi_sxy_y[1][1]), nbyte_pml_y);
            bzero(&(mpw->psi_vyy[1][1]), nbyte_pml_y);
            bzero(&(mpw->psi_vxy[1][1]), nbyte_pml_y);
        }
    }

    /* 4th order */
    if (gv->FDORDER_TIME == 4) {
        bzero(&(mpw->vxx_1[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->vxx_2[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->vxx_3[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->vxx_4[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->vyy_1[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->vyy_2[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->vyy_3[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->vyy_4[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->vxy_1[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->vxy_2[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->vxy_3[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->vxy_4[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->vyx_1[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->vyx_2[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->vyx_3[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->vyx_4[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->svx_1[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->svx_2[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->svx_3[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->svx_4[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->svy_1[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->svy_2[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->svy_3[-gv->ND + 1][-gv->ND + 1]), nbyte_2);
        bzero(&(mpw->svy_4[-gv->ND + 1][-gv->ND + 1]), nbyte_2);

        /* viscoelastic, FDORDER_TIME = 4 */
        if (gv->L) {
            bzero(&(mpw->pr_2[-gv->ND + 1][-gv->ND + 1][1]), nbyte_3);
            bzero(&(mpw->pr_3[-gv->ND + 1][-gv->ND + 1][1]), nbyte_3);
            bzero(&(mpw->pr_4[-gv->ND + 1][-gv->ND + 1][1]), nbyte_3);

            bzero(&(mpw->pp_2[-gv->ND + 1][-gv->ND + 1][1]), nbyte_3);
            bzero(&(mpw->pp_3[-gv->ND + 1][-gv->ND + 1][1]), nbyte_3);
            bzero(&(mpw->pp_4[-gv->ND + 1][-gv->ND + 1][1]), nbyte_3);

            bzero(&(mpw->pq_2[-gv->ND + 1][-gv->ND + 1][1]), nbyte_3);
            bzero(&(mpw->pq_3[-gv->ND + 1][-gv->ND + 1][1]), nbyte_3);
            bzero(&(mpw->pq_4[-gv->ND + 1][-gv->ND + 1][1]), nbyte_3);
        }
    }
}
