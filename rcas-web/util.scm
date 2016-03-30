;;; rcas-web - Web interface for RCAS
;;; Copyright © 2016  Ricardo Wurmus <rekado@elephly.net>
;;; Copyright © 2014  David Thompson <davet@gnu.org>
;;;
;;; This program is free software: you can redistribute it and/or
;;; modify it under the terms of the GNU Affero General Public License
;;; as published by the Free Software Foundation, either version 3 of
;;; the License, or (at your option) any later version.
;;;
;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;;; Affero General Public License for more details.
;;;
;;; You should have received a copy of the GNU Affero General Public
;;; License along with this program.  If not, see
;;; <http://www.gnu.org/licenses/>.

(define-module (rcas-web util)
  #:use-module (ice-9 match)
  #:use-module (srfi srfi-1)
  #:use-module (web request)
  #:use-module (web uri)
  #:export (parse-query-string
            request-path-components
            file-extension
            directory?))

(define (parse-query-string query)
  "Parse and decode the URI query string QUERY and return an alist."
  (let lp ((lst (map uri-decode (string-split query (char-set #\& #\=)))))
    (match lst
      ((key value . rest)
       (cons (cons key value) (lp rest)))
      (() '()))))

(define (request-path-components request)
  (split-and-decode-uri-path (uri-path (request-uri request))))

(define (file-extension file-name)
  (last (string-split file-name #\.)))

(define (directory? filename)
  (string=? filename (dirname filename)))
