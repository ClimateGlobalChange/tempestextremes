# !/bin/sh

# ssh -Y cori.nersc.gov            # Cori
ncdiff test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979010100_1979013123.nc  test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979010100_1979013123_p.nc  test/ERA5_AR_test/ERA5_AR/diff_1979010100_1979013123.nc
ncdiff test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979020100_1979022823.nc  test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979020100_1979022823_p.nc  test/ERA5_AR_test/ERA5_AR/diff_1979020100_1979022823.nc
ncdiff test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979030100_1979033123.nc  test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979030100_1979033123_p.nc  test/ERA5_AR_test/ERA5_AR/diff_1979030100_1979033123.nc
ncdiff test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979040100_1979043023.nc  test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979040100_1979043023_p.nc  test/ERA5_AR_test/ERA5_AR/diff_1979040100_1979043023.nc
ncdiff test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979050100_1979053123.nc  test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979050100_1979053123_p.nc  test/ERA5_AR_test/ERA5_AR/diff_1979050100_1979053123.nc
ncdiff test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979060100_1979063023.nc  test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979060100_1979063023_p.nc  test/ERA5_AR_test/ERA5_AR/diff_1979060100_1979063023.nc
ncdiff test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979070100_1979073123.nc  test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979070100_1979073123_p.nc  test/ERA5_AR_test/ERA5_AR/diff_1979070100_1979073123.nc
ncdiff test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979080100_1979083123.nc  test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979080100_1979083123_p.nc  test/ERA5_AR_test/ERA5_AR/diff_1979080100_1979083123.nc
ncdiff test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979090100_1979093023.nc  test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979090100_1979093023_p.nc  test/ERA5_AR_test/ERA5_AR/diff_1979090100_1979093023.nc
ncdiff test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979100100_1979103123.nc  test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979100100_1979103123_p.nc  test/ERA5_AR_test/ERA5_AR/diff_1979100100_1979103123.nc
ncdiff test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979110100_1979113023.nc  test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979110100_1979113023_p.nc  test/ERA5_AR_test/ERA5_AR/diff_1979110100_1979113023.nc
ncdiff test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979120100_1979123123.nc  test/ERA5_AR_test/ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch.ll025sc.1979120100_1979123123_p.nc  test/ERA5_AR_test/ERA5_AR/diff_1979120100_1979123123.nc

ncview test/ERA5_AR_test/ERA5_AR/diff_1979010100_1979013123.nc
ncview test/ERA5_AR_test/ERA5_AR/diff_1979020100_1979022823.nc
ncview test/ERA5_AR_test/ERA5_AR/diff_1979030100_1979033123.nc
ncview test/ERA5_AR_test/ERA5_AR/diff_1979040100_1979043023.nc
ncview test/ERA5_AR_test/ERA5_AR/diff_1979050100_1979053123.nc
ncview test/ERA5_AR_test/ERA5_AR/diff_1979060100_1979063023.nc
ncview test/ERA5_AR_test/ERA5_AR/diff_1979070100_1979073123.nc
ncview test/ERA5_AR_test/ERA5_AR/diff_1979080100_1979083123.nc
ncview test/ERA5_AR_test/ERA5_AR/diff_1979090100_1979093023.nc
ncview test/ERA5_AR_test/ERA5_AR/diff_1979100100_1979103123.nc
ncview test/ERA5_AR_test/ERA5_AR/diff_1979110100_1979113023.nc
ncview test/ERA5_AR_test/ERA5_AR/diff_1979120100_1979123123.nc

