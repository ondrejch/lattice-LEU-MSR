#!/usr/bin/python3
#
# Test module
#
# Ondrej Chvala, ochvala@utk.edu
# 2019-08-15
# GNU/GPL

import unittest
import salts

my_accuracy = 12    # Digits of accuracy for float tests

class Tests(unittest.TestCase):
    'Class for unit tests'
    def test_flibe_mat(self):
        expected_mat='''\
% Fuel salt: 72%LiF + 16%BeF2 + 12%UF4, U enrichment 0.02
mat fuelsalt  -3.34599934 rgb 240 30 30 tmp  900.000
  3006.09c  -0.000003387223    %  Li-6
  3007.09c  -0.079012792186    %  Li-7
  4009.09c  -0.022555203062    %  Be-9
  9019.09c  -0.451707694555    %  F-19
 92234.09c  -0.000078196964    %  U-234
 92235.09c  -0.008823828716    %  U-235
 92236.09c  -0.000040762584    %  U-236
 92238.09c  -0.437778134710    %  U-238
'''     
        s = salts.Salt()
        self.assertEqual(s.serpent_mat(), expected_mat)

    def test_flibe_rho_700C(self):
        expected_rho = 3.296680695479546
        s = salts.Salt()
        self.assertAlmostEqual(s.densityC(700), expected_rho, my_accuracy)

    def test_flibe_rho_900K(self):
        expected_rho = 3.3459993362645335
        s = salts.Salt()
        self.assertAlmostEqual(s.densityK(900), expected_rho, my_accuracy)

    def test_fnak_mat(self):
        expected_mat='''\
% Fuel salt: 50.5%NaF + 21.5%KF + 28%UF4, U enrichment 0.02
mat fuelsalt  -4.21597326 rgb 240 30 30 tmp  900.000
  9019.09c  -0.287451821554    %  F-19
 11023.09c  -0.095467629714    %  Na-23
 19039.09c  -0.064241464428    %  K-39
 19040.09c  -0.000008266533    %  K-40
 19041.09c  -0.004873891758    %  K-41
 92234.09c  -0.000095917979    %  U-234
 92235.09c  -0.010823486903    %  U-235
 92236.09c  -0.000050000211    %  U-236
 92238.09c  -0.536987520920    %  U-238
'''     
        s = salts.Salt("50.5%NaF + 21.5%KF + 28%UF4")
        self.assertEqual(s.serpent_mat(), expected_mat)

    def test_fnak_rho_700C(self):
        expected_rho = 4.14876616927552
        s = salts.Salt("50.5%NaF + 21.5%KF + 28%UF4")
        self.assertAlmostEqual(s.densityC(700), expected_rho, my_accuracy)

    def test_fnak_rho_900K(self):
        expected_rho = 4.2159732563883845
        s = salts.Salt("50.5%NaF + 21.5%KF + 28%UF4")
        self.assertAlmostEqual(s.densityK(900), expected_rho, my_accuracy)


if __name__ == '__main__':
    unittest.main()        
