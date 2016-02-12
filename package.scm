;;; RCAS -- standalone RNA centric annotation system
;;; Copyright © 2016 Ricardo Wurmus <rekado@elephly.net>
;;;
;;; This file is part of RCAS.
;;;
;;; RCAS is free software; see LICENSE file for details.

;;; Note that RCAS currently depends on the MEME suite, which is
;;; released under a non-free license.  A package definition for MEME
;;; suite is included in this file.  Run the following command to
;;; enter a development environment for RCAS:
;;;
;;;  $ guix environment -l package.scm

(use-modules ((guix licenses) #:prefix license:)
             (guix packages)
             (guix download)
             (guix git-download)
             (guix utils)
             (guix build-system gnu)
             (gnu packages)
             (gnu packages autotools)
             (gnu packages bioinformatics)
             (gnu packages statistics)
             (gnu packages python)

             ;; for meme
             (gnu packages admin)
             (gnu packages base)
             (gnu packages compression)
             (gnu packages perl)
             (gnu packages xml)
             (gnu packages mpi)
             (gnu packages ghostscript)
             (gnu packages web))

(define* (non-free uri #:optional (comment ""))
  "Return a license that does not fit any of the ones above or a collection of
licenses, not approved as free by the FSF.  More details can be found at URI."
  ((@@ (guix licenses) license) "non-free"
           uri
           comment))

;; The latest release of bedtools (version 2.25.0 as of this writing)
;; segfaults on intersect, so we need a more recent version.  We
;; cannot use the latest commit, because that removes command line
;; options from "fastaFromBed".  Another possible version to use is
;; 2.24.x, which does not exhibit this behaviour.
(define-public bedtools/patched
  (let ((commit "2bbe124cbc6e8df393d2ca21b0b29411aefdcf8f"))
    (package
      (inherit bedtools)
      (name "bedtools")
      (version "2.25.0-50-g2bbe124")
      (source (origin
                (method git-fetch)
                (uri (git-reference
                      (url "https://github.com/arq5x/bedtools2.git")
                      (commit commit)))
                (file-name (string-append name "-" version "-checkout"))
                (sha256
                 (base32
                  "1zhzm4kkpldc8l2fiicgb8aywcbcngy21baxpdi3cdcrjqmcfj17"))))
      (arguments
       (substitute-keyword-arguments (package-arguments bedtools)
         ((#:phases phases)
          `(modify-phases ,phases
             (replace 'check
               (lambda _
                 (with-directory-excursion "test"
                   (zero? (system* "bash" "test.sh"))))))))))))

(define meme
  (package
    (name "meme")
    (version "4.11.0")
    (source (origin
              (method url-fetch)
              (uri (string-append "http://meme-suite.org/meme-software/"
                                  version "/meme_" version ".tar.gz"))
              (sha256
               (base32
                "0zr3gvz4k30ggx0wqrbxdrr446vcc1v8q25xnjyjbvqn90gq9i2x"))))
    (build-system gnu-build-system)
    (arguments
     `(#:phases
       (modify-phases %standard-phases
         (add-after 'unpack 'fix-paths-to-tools
           (lambda _
             (substitute* "src/utils.c"
               (("\"hostname")
                (string-append "\"" (which "hostname"))))
             #t))
         (add-after 'unpack 'remove-unused-tests
           (lambda _
             ;; We don't build the web server stuff, so we don't need
             ;; to run the tests for that either.
             (substitute* "tests/scripts/Makefile.in"
               (("tomtom.test") ""))))
         (add-before 'configure 'check-perl-dependencies
           (lambda _
             (zero? (system* "perl" "./scripts/dependencies.pl"))))
         (add-after 'install 'wrap-perl-scripts
           (lambda* (#:key inputs outputs #:allow-other-keys)
             ;; Make sure perl scripts finds all perl inputs at runtime.
             (let ((out (assoc-ref outputs "out")))
               (for-each (lambda (prog)
                           (wrap-program (string-append out "/bin/" prog)
                             `("PERL5LIB" ":" prefix
                               (,(getenv "PERL5LIB")))))
                         '("ama-qvalues"
                           "beeml2meme"
                           "chen2meme"
                           "dreme_xml_to_html"
                           "dreme_xml_to_txt"
                           "elm2meme"
                           "fasta-center"
                           "fasta-fetch"
                           "fasta-grep"
                           "fasta-make-index"
                           "fasta-most"
                           "fasta-subsample"
                           "fasta-unique-names"
                           "hart2meme-bkg"
                           "hartemink2psp"
                           "iupac2meme"
                           "jaspar2meme"
                           "mast_xml_to_html"
                           "mast_xml_to_txt"
                           "matrix2meme"
                           "meme-chip"
                           "meme-rename"
                           "meme_xml_to_html"
                           "nmica2meme"
                           "priority2meme"
                           "psp-gen"
                           "rna2meme"
                           "rsat-retrieve-seq"
                           "rsat-supported-organisms"
                           "scpd2meme"
                           "sites2meme"
                           "taipale2meme"
                           "tamo2meme"
                           "tomtom_xml_to_html"
                           "transfac2meme"
                           "uniprobe2meme"))
              #t))))))
    (inputs
     `(("perl" ,perl)
       ("perl-html-parser" ,perl-html-parser)
       ("perl-html-template" ,perl-html-template)
       ("perl-xml-simple" ,perl-xml-simple)
       ("perl-xml-compile" ,perl-xml-compile)
       ("perl-xml-compile-wsdl11" ,perl-xml-compile-wsdl11)
       ("perl-xml-parser" ,perl-xml-parser)
       ("python" ,python-2) ;only works with Python 2
       ("libxml2" ,libxml2)
       ("libxslt" ,libxslt)
       ("openmpi" ,openmpi)
       ("ghostscript" ,ghostscript)
       ("inetutils" ,inetutils) ;for "hostname"
       ("zlib" ,zlib)))
    (propagated-inputs
     ;; "which" must be propagated because of the weird way it is used
     ;; in "src/exec_parallel.c".  The buffer "cmd_len" is arranged to
     ;; be 6 characters longer than the argument, just enough for the
     ;; string "which ".  I don't want to mess with pointers and
     ;; buffer lengths just to hardcode a path to the "which"
     ;; executable.
     `(("which" ,which)))
    (home-page "http://www.tbi.univie.ac.at/RNA/index.html")
    (synopsis "Motif-based sequence analysis tools")
    (description
     "The MEME Suite allows the biologist to discover novel motifs in
collections of unaligned nucleotide or protein sequences, and to
perform a wide variety of other motif-based analyses.

The MEME Suite supports motif-based analysis of DNA, RNA and protein
sequences.  It provides motif discovery algorithms using both
probabilistic and discrete models, which have complementary strengths.
It also allows discovery of motifs with arbitrary insertions and
deletions (GLAM2).  In addition to motif discovery, the MEME Suite
provides tools for scanning sequences for matches to motifs (FIMO,
MAST and GLAM2Scan), scanning for clusters of motifs (MCAST),
comparing motifs to known motifs (Tomtom), finding preferred spacings
between motifs (SpaMo), predicting the biological roles of
motifs (GOMo), measuring the positional enrichment of sequences for
known motifs (CentriMo), and analyzing ChIP-seq and other large
datasets (MEME-ChIP).")
    (license (non-free "http://meme-suite.org/doc/copyright.html"
                       "license forbids commercial usage"))))

(package
  (name "rcas")
  (version "0.0.0")
  (source #f)
  (build-system gnu-build-system)
  (inputs
   `(("snakemake" ,snakemake)
     ("bedtools" ,bedtools/patched)
     ("meme" ,meme)
     ("python" ,python)
     ("r" ,r)
     ("r-data-table" ,r-data-table)
     ("r-biomart" ,r-biomart)
     ("r-org-hs-eg-db" ,r-org-hs-eg-db)
     ("r-org-ce-eg-db" ,r-org-ce-eg-db)
     ("r-org-dm-eg-db" ,r-org-dm-eg-db)
     ("r-org-mm-eg-db" ,r-org-mm-eg-db)
     ("r-topgo" ,r-topgo)
     ("r-dt" ,r-dt)
     ("r-plotly" ,r-plotly)
     ("r-dplyr" ,r-dplyr)
     ("r-genomation" ,r-genomation)
     ("r-genomicfeatures" ,r-genomicfeatures)
     ("r-rtracklayer" ,r-rtracklayer)
     ("r-rmarkdown" ,r-rmarkdown)))
  (native-inputs
   `(("autoconf" ,autoconf)
     ("automake" ,automake)))
  (synopsis "RNA-centric annotation system")
  (description
   "RCAS aims to be a standalone RNA-centric annotation system that
provides intuitive reports and publication-ready graphics.")
  (home-page "TODO")
  (license license:expat))
