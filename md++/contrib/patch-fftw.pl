#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $fftw_dir = 'fftw-3.3alpha1';
unless(GetOptions('fftw=s' => \$fftw_dir)) {
  print STDERR "Usage: perl $0 --fftw <fftw3 directory to patch>\n";
  exit(1);
}

sub patch_file {
  my $file = shift;
  my $search = shift;
  my $replace = shift;

  # read the file to mem
  open(FILE, $file) or die "Cannot read file $file: $!\n";
  my @lines = <FILE>;
  close(FILE);

  open(FILE, "> $file") or die "Cannot write file $file: $!\n";
  for(my $i = 0; $i <= $#lines; $i++) {
    $lines[$i] =~ s/\Q$search\E/$replace/g;
    print FILE $lines[$i];
  }
  close(FILE);
}

patch_file($fftw_dir . '/api/fftw3.h', 'CONCAT(fftwf_,', 'CONCAT(fftw3f_,');
patch_file($fftw_dir . '/api/fftw3.h', 'CONCAT(fftwl_,', 'CONCAT(fftw3l_,');
patch_file($fftw_dir . '/api/fftw3.h', 'CONCAT(fftw_,', 'CONCAT(fftw3_,');
patch_file($fftw_dir . '/api/fftw3.h', 'FFTW_DEFINE_API(FFTW_MANGLE_DOUBLE, double, fftw_complex)', 'FFTW_DEFINE_API(FFTW_MANGLE_DOUBLE, double, fftw3_complex)');
patch_file($fftw_dir . '/api/fftw3.h', 'FFTW_DEFINE_API(FFTW_MANGLE_FLOAT, float, fftwf_complex)', 'FFTW_DEFINE_API(FFTW_MANGLE_FLOAT, float, fftw3f_complex)');
patch_file($fftw_dir . '/api/fftw3.h', 'FFTW_DEFINE_API(FFTW_MANGLE_LONG_DOUBLE, long double, fftwl_complex)', 'FFTW_DEFINE_API(FFTW_MANGLE_LONG_DOUBLE, long double, fftw3l_complex)');
patch_file($fftw_dir . '/kernel/ifftw.h', 'CONCAT(fftwf_,', 'CONCAT(fftw3f_,');
patch_file($fftw_dir . '/kernel/ifftw.h', 'CONCAT(fftwl_,', 'CONCAT(fftw3l_,');
patch_file($fftw_dir . '/kernel/ifftw.h', 'CONCAT(fftw_,', 'CONCAT(fftw3_,');
patch_file($fftw_dir . '/mpi/fftw3-mpi.h', 'FFTW_MPI_DEFINE_API(FFTW_MPI_MANGLE_DOUBLE, FFTW_MANGLE_DOUBLE, double, fftw_complex)', 'FFTW_MPI_DEFINE_API(FFTW_MPI_MANGLE_DOUBLE, FFTW_MANGLE_DOUBLE, double, fftw3_complex)');
patch_file($fftw_dir . '/mpi/fftw3-mpi.h', 'FFTW_MPI_DEFINE_API(FFTW_MPI_MANGLE_FLOAT, FFTW_MANGLE_FLOAT, float, fftwf_complex)', 'FFTW_MPI_DEFINE_API(FFTW_MPI_MANGLE_FLOAT, FFTW_MANGLE_FLOAT, float, fftw3f_complex)');
patch_file($fftw_dir . '/mpi/fftw3-mpi.h', 'FFTW_MPI_DEFINE_API(FFTW_MPI_MANGLE_LONG_DOUBLE, FFTW_MANGLE_LONG_DOUBLE, long double, fftwl_complex)', 'FFTW_MPI_DEFINE_API(FFTW_MPI_MANGLE_LONG_DOUBLE, FFTW_MANGLE_LONG_DOUBLE, long double, fftw3l_complex)');
patch_file($fftw_dir . '/tests/fftw-bench.h', 'CONCAT(fftwf_,', 'CONCAT(fftw3f_,');
patch_file($fftw_dir . '/tests/fftw-bench.h', 'CONCAT(fftwl_,', 'CONCAT(fftw3l_,');
patch_file($fftw_dir . '/tests/fftw-bench.h', 'CONCAT(fftw_,', 'CONCAT(fftw3_,');
patch_file($fftw_dir . '/tools/fftw-wisdom.c', 'CONCAT(fftwf_,', 'CONCAT(fftw3f_,');
patch_file($fftw_dir . '/tools/fftw-wisdom.c', 'CONCAT(fftwl_,', 'CONCAT(fftw3l_,');
patch_file($fftw_dir . '/tools/fftw-wisdom.c', 'CONCAT(fftw_,', 'CONCAT(fftw3_,');
