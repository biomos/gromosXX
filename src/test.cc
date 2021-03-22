inline void interaction::Perturbed_Nonbonded_Term::
set_lambda_2d( double const lamI_lj, double const lamJ_lj,
            double const lamI_ljs, double const lamJ_ljs,
            double const lamI_crf, double const lamJ_crf,
            double const lamI_crfs, double const lamJ_crfs,
            double const lamI_lj_deriv, double const lamJ_lj_deriv,
            double const lamI_ljs_deriv, double const lamJ_ljs_deriv,
            double const lamI_crf_deriv, double const lamJ_crf_deriv,
            double const lamI_crfs_deriv, double const lamJ_crfs_deriv) {

  // LJ
  lamA_lj = (1-lamI_lj) * (1-lamJ_lj);
  lamAB_lj = (1-lamI_lj) * lamJ_lj;
  lamBA_lj = lamI_lj * (1-lamJ_lj);
  lamB_lj = lamI_lj * lamJ_lj;

  lamA_ljs = (lamI_ljs + lamJ_ljs) * (lamI_ljs + lamJ_ljs);
  lamAB_ljs = (lamI_ljs + (1-lamJ_ljs)) * (lamI_ljs + (1-lamJ_ljs));
  lamBA_ljs = ((1-lamI_ljs) + lamJ_ljs) * ((1-lamI_ljs) + lamJ_ljs);
  lamB_ljs = ((1-lamI_ljs) + (1-lamJ_ljs)) * ((1-lamI_ljs) + (1-lamJ_ljs));

  lamA_lj_deriv = -(1-lamJ_lj)*lamI_lj_deriv - (1-lamI_lj)*lamJ_lj_deriv;
  lamAB_lj_deriv = -lamJ_lj*lamI_lj_deriv + (1-lamI_lj)*lamJ_lj_deriv;
  lamBA_lj_deriv = (1-lamJ_lj)*lamI_lj_deriv - lamI_lj* lamJ_lj_deriv;
  lamB_lj_deriv = lamJ_lj*lamI_lj_deriv + lamI_lj*lamJ_lj_deriv;

  double dfA_dlamI_ljs = 2*(lamI_ljs+lamJ_ljs);
  double dfAB_dlamI_ljs = 2*(1+lamI_ljs-lamJ_ljs);
  double dfBA_dlamI_ljs = -2*(1-lamI_ljs+lamJ_ljs);
  double dfB_dlamI_ljs = -2*(2-lamI_ljs-lamJ_ljs);
  
  double dfA_dlamJ_ljs = dfA_dlamI_ljs;
  double dfAB_dlamJ_ljs = -dfAB_dlamI_ljs;
  double dfBA_dlamJ_ljs = -dfBA_dlamI_ljs;
  double dfB_dlamJ_ljs = dfB_dlamI_ljs;

  lamA_ljs_deriv = dfA_dlamI_ljs * lamI_ljs_deriv + dfA_lamJ_ljs * lamJ_ljs_deriv;
  lamAB_ljs_deriv = dfAB_dlamI_ljs * lamI_ljs_deriv + dfAB_lamJ_ljs * lamJ_ljs_deriv;
  lamBA_ljs_deriv = dfBA_dlamI_ljs * lamI_ljs_deriv + dfBA_lamJ_ljs * lamJ_ljs_deriv;
  lamB_ljs_deriv = dfB_dlamI_ljs * lamI_ljs_deriv + dfB_lamJ_ljs * lamJ_ljs_deriv;

  // CRF
  lamA_crf = (1-lamI_crf) * (1-lamJ_crf);
  lamAB_crf = (1-lamI_crf) * lamJ_crf;
  lamBA_crf = lamI_crf * (1-lamJ_crf);
  lamB_crf = lamI_crf * lamJ_crf;

  lamA_crfs = (lamI_crfs + lamJ_crfs) * (lamI_crfs + lamJ_crfs);
  lamAB_crfs = (lamI_crfs + (1-lamJ_crfs)) * (lamI_crfs + (1-lamJ_crfs));
  lamBA_crfs = ((1-lamI_crfs) + lamJ_crfs) * ((1-lamI_crfs) + lamJ_crfs);
  lamB_crfs = ((1-lamI_crfs) + (1-lamJ_crfs)) * ((1-lamI_crfs) + (1-lamJ_crfs));

  lamA_crf_deriv = -(1-lamJ_crf)*lamI_crf_deriv - (1-lamI_crf)*lamJ_crf_deriv;
  lamAB_crf_deriv = -lamJ_crf*lamI_crf_deriv + (1-lamI_crf)*lamJ_crf_deriv;
  lamBA_crf_deriv = (1-lamJ_crf)*lamI_crf_deriv - lamI_crf* lamJ_crf_deriv;
  lamB_crf_deriv = lamJ_crf*lamI_crf_deriv + lamI_crf*lamJ_crf_deriv;

  double dfA_dlamI_crfs = 2*(lamI_crfs+lamJ_crfs);
  double dfAB_dlamI_crfs = 2*(1+lamI_crfs-lamJ_crfs);
  double dfBA_dlamI_crfs = -2*(1-lamI_crfs+lamJ_crfs);
  double dfB_dlamI_crfs = -2*(2-lamI_crfs-lamJ_crfs);

  double dfA_dlamJ_crfs = dfA_dlamI_crfs;
  double dfAB_dlamJ_crfs = -dfAB_dlamI_crfs;
  double dfBA_dlamJ_crfs = -dfBA_dlamI_crfs;
  double dfB_dlamJ_crfs = dfB_dlamI_crfs;

  lamA_crfs_deriv = dfA_dlamI_crfs * lamI_crfs_deriv + dfA_lamJ_crfs * lamJ_crfs_deriv;
  lamAB_crfs_deriv = dfAB_dlamI_crfs * lamI_crfs_deriv + dfAB_lamJ_crfs * lamJ_crfs_deriv;
  lamBA_crfs_deriv = dfBA_dlamI_crfs * lamI_crfs_deriv + dfBA_lamJ_crfs * lamJ_crfs_deriv;
  lamB_crfs_deriv = dfB_dlamI_crfs * lamI_crfs_deriv + dfB_lamJ_crfs * lamJ_crfs_deriv;

}
          


inline void interaction::Perturbed_Nonbonded_Term
::lj_crf_soft_2d_interaction(math::Vec const &r,
        double const A_i_c6, double const A_i_c12,
        double const B_i_c6, double const B_i_c12,
        double const A_j_c6, double const A_j_c12,
        double const B_j_c6, double const B_j_c12,
        double const A_i_q, double const B_i_q,
        double const A_j_q, double const B_j_q,
        double const alpha_lj, double const alpha_crf,
        double &force1, double &force6, double &force12,
        double &e_lj, double &e_crf,
        double &de_lj, double & de_crf, unsigned int eps) {


  //LJ
  const double A_c6 = A_i_c6 * A_j_c6;
  const double A_c12 = A_i_c12 * A_j_c12;

  const double AB_c6 = A_i_c6 * B_j_c6;
  const double AB_c12 = A_i_c12 * B_j_c12;

  const double BA_c6 = B_i_c6 * A_j_c6;
  const double BA_c12 = B_i_c12 * A_j_c12;

  const double B_c6 = B_i_c6 * A_j_c6;
  const double B_c12 = B_i_c12 * B_j_c12;

  double A_c126, AB_c126, BA_c126, B_c126;
  if (A_c6 != 0) A_c126 = A_c12 / A_c6;
  else A_c126 = 0.0;
  if (AB_c6 != 0) AB_c126 = AB_c12 / AB_c6;
  else AB_c126 = 0.0;
  if (BA_c6 != 0) BA_c126 = BA_c12 / BA_c6;
  else BA_c126 = 0.0;
  if (B_c6 != 0) B_c126 = B_c12 / B_c6;
  else B_c126 = 0.0;

  const double dist4 = dist2 * dist2;
  const double dist6 = dist4 * dist2;

  const double A_dist6soft = dist6 + alpha_lj * lamA_ljs*A_c126;
  const double AB_dist6soft= dist6 + alpha_lj * lamAB_ljs*AB_c126;
  const double BA_dist6soft= dist6 + alpha_lj * lamBA_ljs*BA_c126;
  const double B_dist6soft = dist6 + alpha_lj * lamB_ljs*B_c126;

  const double A_dist6isoft = 1.0 / A_dist6soft;
  const double AB_dist6isoft = 1.0 / A_dist6soft;
  const double BA_dist6isoft = 1.0 / B_dist6soft;
  const double B_dist6isoft = 1.0 / B_dist6soft;

  force6 = -6.0 * (lamA_lj * A_c6 * A_dist6isoft * A_dist6isoft +
          lamAB_lj * AB_c6 * AB_dist6isoft * AB_dist6isoft +
          lamBA_lj * BA_c6 * BA_dist6isoft * BA_dist6isoft +
          lamB_lj * B_c6 * B_dist6isoft * B_dist6isoft) * dist4;

  force12 = 12 * (lamA_lj * A_c12 * A_dist6isoft * A_dist6isoft * A_dist6isoft +
          lamAB_lj * AB_c12 * AB_dist6isoft * AB_dist6isoft * AB_dist6isoft +
          lamBA_lj * BA_c12 * BA_dist6isoft * BA_dist6isoft * BA_dist6isoft +
          lamB_lj * B_c12 * B_dist6isoft * B_dist6isoft * B_dist6isoft) * dist4;

  const double A_e_lj = (A_c12 * A_dist6isoft - A_c6) * A_dist6isoft;
  const double B_e_lj = (B_c12 * B_dist6isoft - B_c6) * B_dist6isoft;
  const double AB_e_lj = (AB_c12 * AB_dist6isoft - AB_c6) * AB_dist6isoft;
  const double BA_e_lj = (BA_c12 * BA_dist6isoft - BA_c6) * BA_dist6isoft;

  e_lj = lamA_lj * A_e_lj + lamAB_lj * AB_e_lj + 
         lamBA_lj * BA_e_lj + lamB_lj * B_e_lj;

  de_lj = A_e_lj * lamA_lj_deriv - (alpha_lj * lamA_lj * A_c126 * A_dist6isoft * A_dist6isoft) *
          (2 * A_c12 * A_dist6isoft - A_c6) * lamA_ljs_deriv +
          AB_e_lj * lamAB_lj_deriv - (alpha_lj * lamAB_lj * AB_c126 * AB_dist6isoft * AB_dist6isoft) *
          (2 * AB_c12 * AB_dist6isoft - AB_c6) * lamAB_ljs_deriv +
          BA_e_lj * lamBA_lj_deriv - (alpha_lj * lamBA_lj * BA_c126 * BA_dist6isoft * BA_dist6isoft) *
          (2 * BA_c12 * BA_dist6isoft - BA_c6) * lamAB_ljs_deriv +
          B_e_lj * lamB_lj_deriv - (alpha_lj * lamB_lj * B_c126 * B_dist6isoft * B_dist6isoft) *
          (2 * B_c12 * B_dist6isoft - B_c6) * lamB_ljs_deriv;


  //CRF
  const double dist2 = abs2(r);
  assert(dist2 != 0);

  double A_q = A_i_q * A_j_q;
  double AB_q = A_i_q * B_j_q;
  double BA_q = B_i_q * A_j_q;
  double B_q = B_i_q * B_j_q;

  const double A_dist2soft = dist2 + alpha_crf*lamA_crfs;
  const double AB_dist2soft = dist2 + alpha_crf*lamAB_crfs;
  const double BA_dist2soft = dist2 + alpha_crf*lamBA_crfs;
  const double B_dist2soft = dist2 + alpha_crf*lamB_crfs;

  const double A_distisoft = 1.0 / sqrt(A_dist2soft);
  const double AB_distisoft = 1.0 / sqrt(AB_dist2soft);
  const double BA_distisoft = 1.0 / sqrt(BA_dist2soft);
  const double B_distisoft = 1.0 / sqrt(B_dist2soft);

  const double A_dist3isoft = A_distisoft / A_dist2soft;
  const double AB_dist3isoft = AB_distisoft / AB_dist2soft;
  const double BA_dist3isoft = BA_distisoft / BA_dist2soft;
  const double B_dist3isoft = B_distisoft / B_dist2soft;

  const double A_cut2soft = m_cut2 + alpha_crf * lamA_crfs;
  const double AB_cut2soft = m_cut2 + alpha_crf * lamAB_crfs;
  const double BA_cut2soft = m_cut2 + alpha_crf * lamBA_crfs;
  const double B_cut2soft = m_cut2 + alpha_crf * lamB_crfs;

  const double A_cut2soft3 = A_cut2soft * A_cut2soft * A_cut2soft;
  const double AB_cut2soft3 = AB_cut2soft * AB_cut2soft * AB_cut2soft;
  const double BA_cut2soft3 = BA_cut2soft * BA_cut2soft * BA_cut2soft;
  const double B_cut2soft3 = B_cut2soft * B_cut2soft * B_cut2soft;

  const double A_crf_2cut3i = m_crf_2[eps] / sqrt(A_cut2soft3);
  const double AB_crf_2cut3i = m_crf_2[eps] / sqrt(AB_cut2soft3);
  const double BA_crf_2cut3i = m_crf_2[eps] / sqrt(BA_cut2soft3);
  const double B_crf_2cut3i = m_crf_2[eps] / sqrt(B_cut2soft3);

  const double A_crf_cut3i = 2 * A_crf_2cut3i;
  const double AB_crf_cut3i = 2 * AB_crf_2cut3i;
  const double BA_crf_cut3i = 2 * BA_crf_2cut3i;
  const double B_crf_cut3i = 2 * B_crf_2cut3i;

  const double A_crf_pert = 3.0 * A_crf_2cut3i / A_cut2soft;
  const double AB_crf_pert = 3.0 * AB_crf_2cut3i / AB_cut2soft;
  const double BA_crf_pert = 3.0 * BA_crf_2cut3i / BA_cut2soft;
  const double B_crf_pert = 3.0 * B_crf_2cut3i / B_cut2soft;

  force1 = (lamA_crf * A_q * (A_dist3isoft + A_crf_cut3i) +
          lamAB_crf * AB_q * (AB_dist3isoft + AB_crf_cut3i) +
          lamBA_crf * BA_q * (BA_dist3isoft + BA_crf_cut3i) +
          lamB_crf * B_q * (B_dist3isoft + B_crf_cut3i)) *
          math::four_pi_eps_i;

  const double A_e_crf = A_q * (A_distisoft - A_crf_2cut3i * dist2 - m_crf_cut[eps]);
  const double AB_e_crf = AB_q * (AB_distisoft - AB_crf_2cut3i * dist2 - m_crf_cut[eps]);
  const double BA_e_crf = BA_q * (BA_distisoft - BA_crf_2cut3i * dist2 - m_crf_cut[eps]);
  const double B_e_crf = B_q * (B_distisoft - B_crf_2cut3i * dist2 - m_crf_cut[eps]);

  e_crf = (lamA_crf * A_e_crf + lamAB_crf * AB_e_crf + 
           lamBA_crf * BA_e_crf + lamB_crf * B_e_crf) * math::four_pi_eps_i;

  de_crf = ((A_e_crf * lamA_crf_deriv - (alpha_crf /2) * lamA_crf *
           (A_dist3isoft - A_crf_pert * dist2) * lamA_crfs_deriv ) * A_q +
           (AB_e_crf * lamAB_crf_deriv - (alpha_crf /2) * lamAB_crf *
           (AB_dist3isoft - AB_crf_pert * dist2) * lamAB_crfs_deriv ) * AB_q +
           (BA_e_crf * lamBA_crf_deriv - (alpha_crf /2) * lamBA_crf *
           (BA_dist3isoft - BA_crf_pert * dist2) * lamBA_crfs_deriv ) * BA_q  +
           (B_e_crf * lamB_crf_deriv - (alpha_crf /2) * lamB_crf *
           (B_dist3isoft - B_crf_pert * dist2) * lamB_crfs_deriv ) * B_q ) 
           * math::four_pi_eps_i;

}

