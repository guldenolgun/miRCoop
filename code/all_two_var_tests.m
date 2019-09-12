%Dino Sejdinovic, 2013
%D. Sejdinovic, A. Gretton and W. Bergsma.  A KERNEL TEST FOR THREE-VARIABLE INTERACTIONS, 2013.
%runs all two-variable tests for a three-variable dataset
function TwoVar = all_two_var_tests(xx,yy,zz,param)
    [TwoVar.XY.pvalue, TwoVar.XY.rejectNull]     = two_var_test(xx,yy,param);
    [TwoVar.XZ.pvalue, TwoVar.XZ.rejectNull]     = two_var_test(xx,zz,param);
    [TwoVar.YZ.pvalue, TwoVar.YZ.rejectNull]     = two_var_test(yy,zz,param);
    [TwoVar.XY_Z.pvalue, TwoVar.XY_Z.rejectNull] = two_var_test([xx yy],zz,param);
    [TwoVar.XZ_Y.pvalue, TwoVar.XZ_Y.rejectNull] = two_var_test([xx zz],yy,param);
    [TwoVar.YZ_X.pvalue, TwoVar.YZ_X.rejectNull] = two_var_test([yy zz],xx,param);
end

