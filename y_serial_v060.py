#  y_serial Python module         Version 0.60                Date : 2010-08-20
# -*- coding: iso-8859-1 -*-
#                                 http://yserial.sourceforge.net
''' 
_______________ y_serial :: warehouse compressed Python objects with SQLite

  Dependencies:  at least Python v2.5 because it includes the sqlite3 module. 
                 Database itself and all other modules used are standard issue.


          _____ PREFACE (in reStructured Text format for our site)

Intro and quick example
-----------------------

*In about ten minutes, you should be able to simply give some
labels to any Python object and save it to a database file; then
get back a set of objects by specifying portion of their labels.*

Here's a quick EXAMPLE for a typical situation. After downloading
the y_serial module,  create an instance which is associated with a
regular file:: 

     import y_serial_v060 as y_serial
     demo = y_serial.Main( '/tmp/agency.sqlite' )
     #  ^instance          ^include suitable path for database file

     #         ... now we do some work producing an object obj, e.g.  
     obj = 911

That object could have been a dictionary with a complex structure,
but let's continue on for the insert::

     demo.insert( obj, "#plan agent007 #london", 'goldfinger' )
     #                 ^notes                    ^table

We label each object with "notes" which can be arbitrarily long
text (or UTF-8), containing keywords or tags, but excluding commas.
Within that file we specify a "table" (merely an organizational 
subsector).  Some time later, perhaps in another script, we will
want to select some object:: 

     eg1 = demo.select( "agent00[1-7],#plan", 'goldfinger' )
     #                   ^search values are space-sensitive
     #                                  and comma separated;
     #                    arbitrarily many permitted in string.

     print "Example 1:  ", eg1
     #            That reveals the _latest_ goldfinger plan 
     #            which involves any one of the top seven agents
     #            anywhere in the world including London.


That's it... **only a few lines of Python code to store compressed
serialized objects in a database, and to selectively retrieve
them** (with optional regular expression, and remarkably, *without
writing any SQL commands*). DEAD SIMPLE -- *with only one module
imported*. Hopefully you see how widely this is applicable...


Installation and license 
------------------------

Sole requirement: Python version 2.x where x is 5 or greater.

Download the latest version of the module at
`http://sourceforge.net/projects/yserial
<http://sourceforge.net/projects/yserial/>`_ and put it where your
Python can find it. No tar.gz or eggs here ;-) The module includes
the tutorial documentation within itself. You are free to use
*y_serial* under the BSD license. No monetary charge at all;
however, if you are a developer, please critically review the code.
More eyeballs lead to increased scrutiny, and thus greater
reliability.

     
Overview 
--------

The purpose of y_serial is to keep data persistent. It is based on
key/value where the conceptual key is 

     - filename + table_name + primary_key + timestamp + notes 

and the value is some object. (Ironically this is implemented using
a SQL database. ;-)  **Our goal is to have the simplest possible
user interface in the foreground, yet have near optimal computing
in the background.**  

By "objects" we generally mean Python objects, but we include
support for files (binary, image, etc.) and URL content (e.g.
webpages).  Python objects are strings, dictionaries, lists,
tuples, classes, and instances.  Objects are inserted into a
hierarchy: database file, table within that database, row within
table.  Moreover, each object is annotated for reference by "notes"
and automatically timestamped.

*So what is happening in the background?* To minimize storage size,
each object is compressed. But before compression, we serialize the
object. The processing between y_serial and the database is handled
by the sqlite3 module **(all dependencies are standard issue
modules)**.  Your program's interaction with y_serial normally will
not require writing SQL, although subqueries can be customized. 

y_serial is written as a single Python module which reads like a
working tutorial and includes many tips and references.  It's
instructive in the way it unifies the standard batteries:

     - sqlite3 (as of Python v2.5)
     - zlib (for compression)
     - cPickle (for serializing objects)
     - ["re" module is not used, instead we access much faster 
       SQLite functions for regular expressions]

Technicalities aside, you are spared from explicitly spelling out
many of the difficult protocol details: cursor/connection,
SQL/DB-API implementation, serialization, compression, search
algorithm, etc. -- for these are optimized to interact with speed,
security, and concurrency -- yet handled transparently. 

Our module is faster than comparable approaches under PostgreSQL.
Among serialization methods, we found cPickle to be one of the
fastest, and so we have used it in a secure manner. Try y_serial
with a few million objects. 

We recommend SQLite because it requires neither separate
installation nor a server process; also, it uses single normal
files (easy to backup or send), not an elaborate filesystem.
Moreover, in comparison to similar applications with MySQL or
PostgreSQL, SQLite is extremely fast and suits most purposes
wonderfully. Should you later decide to migrate out of SQLite,
y_serial can help port your objects elsewhere including other NoSQL
implementations.

*The means for insertion, organization by annotation, and finally
retrieval are designed to be* **simple to use**.  The notes
effectively label the objects placed into the database.  We can
then later query the database, for example, by *regex* (regular
expression) searching on notes, and placing the qualified objects
in a dictionary.  The keys of this dictionary correspond to the
unique primary keys used in the database. If necessary we can
access the timestamp for each object.  We can thus use Python code
to process the contents of the qualified dictionary, in effect, a
data subset.  If the objects in that dictionary are themselves
dictionaries we are essentially dealing with *schema-less data*
(see the compelling Friendfeed case study in the module's
Endnotes).

To illustrate, let's continue our example by adding an object for
agent006::

     obj = 411
     demo.insert( obj, "agent006 #paris #plan", 'goldfinger' )
     
     #     We now can get a dictionary of objects
     #     which matches our search values:
     #
     eg2 = demo.selectdic( "agent00[1-7],#plan", 'goldfinger' )
     print "Example 2:  ", eg2
     #
     #          which should look like:
     #     {1: [1257874696, u'#plan agent007 #london', 411], 
     #      2: [1257874696, u'agent006 #paris #plan', 911]   }
     
Notice that we used a different method called *selectdic* which
produces a dictionary whose keys are the unique primary keys
automatically assigned in the database.  Inside the list are the
(unix) epoch timestamp, followed by (unicode) notes, then object.
This means that we can *work with flexible data subsets using
Python code rather than cumbersome SQL.* 
    

Other features
--------------

Instead of using comma-separated values, as in our example so far,
we could have crafted a custom subquery and used a method called
*dicsub*.

Or we could just skip any subquery altogether. Here we pick out the
most recent n-th entry::

     eg3 = demo.select( 0, 'goldfinger' )
     print "Example 3:  ", eg3

The method called "*view*" will verbosely pretty-print deeply
nested structures::

     demo.view( 5, 'goldfinger' )
     #          ^last m inserts (or use search string argument).

y_serial can also act like a persistent *QUEUE*.  Whatever that is
retrieved can be deleted thereafter by appending "POP=True" at the
end of any applicable method::

     eg4 = demo.select( 0, 'goldfinger', POP=True )

Object(s) can of course be deleted directly::

     demo.delete( "agent00?", 'goldfinger' )
     #                    ^where notes mention any single digit agent.

To get rid of stale data we could freshen a table and vacuum the
entire database via *clean*::

     demo.clean( 365.5 , 'goldfinger' )
     #           ^retain one year-old or less prior to last insert.

To delete the entire table::

     demo.droptable( 'goldfinger' )

Other useful methods are available:

     - insert any external file (via *infile*). This is handy for
       working with thousands of image files.

     - insert anything on the web by URL (via *inweb*).

     - insert in batches by generator (via *ingenerator*). This can
       be used to rapidly store a series of computationally intense 
       results for quick retrieval at a later time.  

For concurrency we can easily code for a farm of databases using
the module's copy functions. In general, your program can have
multiple interacting instances which control distinct database
files.

[*What's in beta?*  In heavy concurrent situations (say, hundreds
of near-simultaneous writes per second), SQLite has a problem
because of the way it locks the database.  We can alleviate the
problem by diffusing such operations across many databases (called
"barns" via class *Farm*), and then writing back to the target
database as a batch (which is far more efficient than single
writes) over time.  The harvest time to reap the batch is
stochastic (see the "plant" method). That introduces some latency
in accessing the newest objects -- the cost for scaling up
concurrency.]

The class *Main* should be stable for all practical purposes.
If you run across any bugs, please kindly report them to the
`Tracker
<https://sourceforge.net/tracker/?group_id=277002&atid=1176411>`_.
For group discussions, check out the SourceForge link in the
left sidebar.

For other specialized features, please RTM "read the module" for
tips on usage.


Summary
-------

**y_serial = serialization + persistance. In a few lines of code,
compress and annotate Python objects into SQLite; then later
retrieve them chronologically by keywords without any SQL. Highly
useful NoSQL "standard" module for a database to store schema-less
data.**


_______________ TABLE OF CONTENTS:

   - Preface with quick example
   - Usage with SUMMARY of CLASSES, METHODS, and FUNCTIONS:

     Base:
          _______________ Attributes and methods for database setup.
               Set path to database for all instances; db0 is default.
               Connection and execution methods.
     Insertion( Base ):
          _______________ INSERT pz BLOB into DATABASE
          inbatch( self, objseq, table=Base.tab0 ):
               Pickle and compress sequence of annotated objects; insert.
       *  ingenerator( self, generate_objnotes, table=Base.tab0 ):
               Pickle and compress via generator function, then insert.
      **  insert( self, obj, notes='#0notes', table=Base.tab0 ):
               Pickle and compress single object; insert with annotation.
     Annex( Insertion ):
          _______________ Add macro-objects (files, URL content) to DATABASE
          inweb( self, URL, notes='', table=Base.tab0 ):
               Pickle and compress URL content, then insert into table.
       *  infile( self, filename, notes='', table=Base.tab0 ):
               Pickle and compress any file, then insert contents into table.
     Answer( Base ):
          _______________ Single item answer shouted out.
          shout( self, question, table=Base.tab0 ):
               Shout a question; get a short answer.
      **  lastkid( self, table=Base.tab0 ):
               Get primary key ID of the last insert.
          lastsec( self, table=Base.tab0 ):
               Get time in unix seconds of the last insert.
          lastdate( self, table=Base.tab0 ):
               Get local date/time of the last insert.
     Util:
          _______________ Utility methods for keys, subquery, comma
          comma2list( self, csvstr, wild=True ):
               Convert comma separated values to a parameter list.
     Deletion( Base, Util ):
          _______________ Deletion methods; also used for queue POP 
          deletesub( self, subquery, parlist=[], table=Base.tab0 ):
               Delete row(s) matching the subquery.
          deletekid( self, kid, table=Base.tab0 ):
               Delete single row with primary key kid.
          deletecomma( self, csvstr, table=Base.tab0, wild=True ):
               Delete row(s): notes match comma separated values in string.
       *  delete( self, dual, table=Base.tab0, wild=True ):
               Alias "delete":  deletekid OR deletecomma.
          droptable( self, table=Base.tab0 ):
               Delete a table: destroys its structure, indexes, data.
     Subquery( Util, Answer, Deletion ):
          _______________ SUBQUERY table, get dictionary. POP QUEUE.
          dicsub(self, subquery='', parlist=[], table=Base.tab0, POP=False):
               Subquery table to get objects into response dictionary.
          diclast( self, m=1, table=Base.tab0, POP=False ):
               Get dictionary with last m consecutive kids in table.
          diccomma( self, csvstr, table=Base.tab0, wild=True, POP=False ):
               Get dictionary where notes match comma separated values.
       *  selectdic( self, dual=1, table=Base.tab0, POP=False ):
               Alias "selectdic":  diclast  OR diccomma.
     Display( Subquery ):
          _______________ View subquery via pretty print
          viewsub(self, subquery='', parlist=[], table=Base.tab0, POP=False):
               Subquery, order keys, and print qualified dictionary.
          viewlast( self, m=1, table=Base.tab0, POP=False ):
               Print last m consecutive kids in table.
          viewcomma(self, csvstr='', table=Base.tab0, wild=True, POP=False):
               Print dictionary where notes match comma separated values.
       *  view( self, dual=1, table=Base.tab0, POP=False ):
               Alias "view":  viewlast OR viewcomma.
     Latest( Display ):
          _______________ Retrieve the latest qualified object "omax" 
          omaxsub(self, subquery='', parlist=[], table=Base.tab0, POP=False):
               Get the latest object omax which matches subquery.
          omaxlast( self, n=0, table=Base.tab0, POP=False ):
               Most quickly get the latest n-th object using key index.
          omaxcomma( self, csvstr, table=Base.tab0, wild=True, POP=False ):
               Get latest object where notes match comma separated values.
      **  select( self, dual=0, table=Base.tab0, POP=False ):
               Alias "select":  omaxlast OR omaxcomma.
       *  getkid( self, kid, table=Base.tab0, POP=False ):
               Retrieve a row given primary key kid, POP optional.
     Oldest( Latest ):
          _______________ Retrieve the oldest qualified object "omin" 
          ominfirst( self, n=0, table=Base.tab0, POP=False ):
               Most quickly get the oldest n-th object using key index.
       *  fifo( self, table=Base.tab0 ):
               FIFO queue: return oldest object, then POP (delete) it.
     Care( Answer, Deletion ):
          _______________ Maintenance methods
          freshen( self, freshdays=None, table=Base.tab0 ):
               Delete rows in table over freshdays-old since last insert.
          vacuum( self ):
               Defrag entire database, i.e. all tables therein.
                    - why VACUUM?
       *  clean( self, freshdays=None, table=Base.tab0 ):
               Delete stale rows after freshdays; vacuum/defrag database.
     Main( Annex, Oldest, Care ):
          _______________ Summary for use of a single database.
     copysub( subquery, parlist, tablex, tabley, dbx=Base.db0, dby=Base.db0 ):
          Subselect from tablex, then copy to tabley (in another database).
     copylast( m, tablex, tabley, dbx=Base.db0, dby=Base.db0 ):
          Copy last m consecutive kids in tablex over to tabley.
   * comma( *string_args ):
          Join string-type arguments with comma (cf. csvstr, comma2list).
     copycomma( csvstr, tablex, tabley, dbx=Base.db0, dby=Base.db0, wild=True ):
          Subselect by comma separated values from tablex, then copy to tabley.
  ** copy( dual, tablex, tabley, dbx=Base.db0, dby=Base.db0, wild=True ):
          Alias "copy":  copylast OR copycomma

     Farm:
          _______________ Start a farm of databases for concurrency and scale.
          __init__( self, dir=dir0, maxbarns=barns0 ):
               Set directory for the farm of maxbarns database files.
          farmin( self, obj, notes='#0notes', table=Base.tab0, n=0 ):
               Insert an object with notes into barn(n).
          harvest(self, dual, tablex,tabley, n, dby=Base.db0, wild=True, size=10):
               After farmin, reap dual under tablex barn(n) by size expectation.
          cleanfarm( self, freshdays=None, table=Base.tab0 ):
               Delete stale rows after freshdays; vacuum/defrag barns.
        * plant( self, obj, notes='#0notes', table=Base.tab0, dby=Base.db0 ):
               FARM SUMMARY: farmin insert with generic self-cleaning harvest.

   - Change log
   - TODO list
   - ENDNOTES: operational tips and commentary with references
        - pz Functions for FILE.gz
             - Database versus FILE.gz
   - tester( database=Base.db0 ):
          Test class Main for bugs. 
   - testfarm( dir=Farm.dir0, maxbarns=Farm.barns0 ):
          Test class Farm for bugs. Include path for directory.
   - Acknowledgements and Revised BSD LICENCE


_______________  CHANGE LOG

     2010-08-20  v0.60:  Certified to run under Python 2.6 series.
                         Edited the preface (also used for welcome page).
                         Added getkid. Cosmetic touch-up for view method.

     2010-04-25  v0.53:  Be sure to assign default database db0 in class Base.

     2009-11-22  v0.52:  Added plant (insert) method summarizing class Farm.

     2009-11-19  v0.51:  Added fifo (queue) method to class Oldest.
                         Added beta Farm (of databases!) for concurrency.

  !! 2009-11-11  v0.50:  TOTAL TOP-DOWN GENERALIZED CLASS REWRITE
                    of the final bottom-up all-function revision 28.

          - Python objects now are warehoused among various database files 
               by design. Each instance of the class Main is associated with 
               a particular database file.
          - Dropped *glob functions in favor of the *comma form (see 
               comma2list to see how this works exactly). 
          - Simplified the names of many functions since most of them 
               are now methods within some class. 
          - LIFO has been renamed as POP to reflect queue generlization.

     2009-10-26  Revision 28: final bottom-up all-function version 
                    which also archives prior change log.

  !! = indicates change(s) which broke backward compatibility.
  pz = pickled zlib compressed binary format.


_______________ TODO List

The sqlite3 module, maintained by Gerhard Haering, has been updated 
from version 2.3.2 in Python 2.5 to version 2.4.1 in Python 2.6. 

[ ] - update code for threads when the following becomes public:

     Contrary to popular belief, newer versions of sqlite3 do support 
     access from multiple threads.  This can be enabled via optional 
     keyword argument "check_same_thread": 

          sqlite.connect(":memory:", check_same_thread = False)

     However, their docs omit this option.  --Noted 2010-05-24  

[ ] - update code for Python version 3.x

     When we run y_serial.tester() with the "-3" flag under Python 2.6, 
     there is only one message:

          "DeprecationWarning: buffer() not supported in 3.x"

     This actually is about how BLOBs are handled in the sqlite3 module, 
     and is being resolved by Gerhard Haering, 
     see http://bugs.python.org/issue7723
     Thereafter, we expect an easy transition to Python 3.
'''

#  _______________ Variable settings with imports

DEBUG    = False
#  Here's how to EXECUTE TESTS. First, be sure to change the default 
#  database file to suit yourself; see assignment db0 in class Base.
#            import y_serial_v053 as y_serial
#            y_serial.tester( database )
#            #        ^for the principal class Main
#            y_serial.testfarm( directory, maxbarns )
#                     ^for the beta version, not yet in Main.

import sqlite3 as ysql
#                 ^ for portability to another SQL database.

import pprint
#  pretty print to Display nested arbitrary Python data structures.

#       __________ pz FUNCTIONS for pickled compression

import cPickle as yPickle
#      ^ written in C, compatible with pickle but thousand times faster.
#        But some situations may only allow pure-Python pickle module.
#        The data stream produced by pickle and cPickle are identical.
#        So take your pick for yPickle.
#        [Future note:  pickle in Python v3 will integrate cPickle.]

pickle_protocol = 2
#  0 = original ASCII protocol and is backwards compatible.
#  1 = old binary format which is also backwards compatible.
#  2 = introduced in Python 2.3, more efficient pickling of new-style classes.
#  3 = for Python 3 with explicit support for bytes objects and byte arrays,
#         but it is not backward compatible with protocol 2.
#         So use that since Python 3 will still understand what to do.

import zlib
#      ^zlib compression for pickled items.

compress_level  = 7
#               1 to 9 controlling the level of compression; 
#       1 is fastest and produces the least compression, 
#       9 is slowest and produces the greatest compression. 

def pzdumps( obj ):
     '''Pickle object, then compress the pickled.'''
     return zlib.compress( yPickle.dumps(obj, pickle_protocol), compress_level)
     #      as binary string.

def pzloads( pzob ):
     '''Inverse of pzdumps:  decompress pz object, then unpickle.'''
     return yPickle.loads( zlib.decompress( pzob ) )  
     #              ^auto-detects pickle protocol.



class Base:
     '''_______________ Essential attributes and methods for database setup.'''

     db0      = '/home/yaya/var/db/y_serial.sqlite'
     #           ========================================================== ***
     #          ^ be sure to specify an absolute path to the database file.
     #            This is just your convenient DEFAULT DATABASE FILE.
     #            Specify other such files explicitly when creating instances.
     #
     #            [ Using an in-memory database ':memory:' will not work here 
     #              because we go in and out of connection as needed. ]

     tab0     = 'tmptable'
     #          ^default SQL table for storing objects temporarily.

     TRANSACT = 'IMMEDIATE'
     #  among: None (autocommit), 'DEFERRED' (default), 'IMMEDIATE', 'EXCLUSIVE'. 
     #
     #  We want to support transactions among multiple concurrent sessions, and
     #  to AVOID DEADLOCK.  For our purposes, such sessions should start out by:
     #       "BEGIN IMMEDIATE TRANSACTION;"
     #  to guarantee a (reserved) write lock while allowing others to read.
     #        See _The Definitive Guide to SQLite_, chapters 4 and 5, 
     #        by Michael Owens (2006, Apress) for the clearest explanation.

     TIMEOUT  = 14
     #          ^ in seconds (default: 5)
     #  During multiple concurrent sessions, writing creates a certain type 
     #  of lock until a transaction is committed. TIMEOUT specifies how long 
     #  a connection should wait for that lock to be released until raising 
     #  an exception.  Increase the wait if a very large amount of objects 
     #  is routinely inserted during a single session.

     def __init__( self, db=db0 ):
          '''Set path to database for all instances; db0 is default.'''
          self.db = db

     def proceed( self, sql, parlist=[[]] ):
          '''Connect, executemany, commit, then finally close.'''
          try:
               con = ysql.connect( self.db,    timeout = self.TIMEOUT, 
                                       isolation_level = self.TRANSACT )
               cur = con.cursor()
               cur.executemany( sql, parlist )
               #        for an empty ^parameter list, use [[]].
               con.commit()
               #   ^MUST remember to commit! else the data is rolled back! 
          except:
               a = " !! Base.proceed did not commit. [Check db path.] \n"
               b = "             Suspect busy after TIMEOUT,          \n"
               c = "             tried this sql and parameter list:   \n"
               raise IOError, "%s%s%s%s\n%s" % ( a, b, c, sql, parlist )
          finally:
               cur.close()
               con.close()
               #   ^ very important to release lock for concurrency.

     def respond( self, klass, sql, parlist=[] ):
          '''Connect, execute select sql, get response dictionary.'''
          try:
               con = ysql.connect( self.db,    timeout = self.TIMEOUT, 
                                       isolation_level = self.TRANSACT )
               cur = con.cursor()
               response = {}
               for tupler in cur.execute( sql, parlist ):
                    self.responder( klass, tupler, response )
               #         ^ to be defined in a subclass
               #           (mostly to process output from subqueries).
               #  con.commit() intentionally omitted.
          except:
               a = " !! Base.respond choked, probably because     \n"
               b = "             object feels out of context.     \n"
               c = "           Tried this sql and parameter list: \n"
               raise IOError, "%s%s%s%s\n%s" % (a, b, c, sql, parlist)
          finally:
               cur.close()
               con.close()
          return response

     def createtable( self, table=tab0 ):
          '''Columns created: key ID, unix time, notes, and pzblob.'''
          a = 'CREATE TABLE IF NOT EXISTS %s' % table
          b = '(kid INTEGER PRIMARY KEY, tunix INTEGER,'
          c = 'notes TEXT, pzblob BLOB)'
          sql = ' '.join( [a, b, c] )
          try:
               self.proceed( sql )
               #  try construct is useful for portability in standard 
               #  cases where clause "IF NOT EXISTS" is not implemented.
          except IOError:
               if DEBUG:
                    print " :: createtable: table exists."
          #    createtable is designed to be harmless if it 
          #    left sitting in your script.



class Insertion( Base ):
     '''_______________ INSERT pz BLOB into DATABASE'''

     #  For inbatch we shall assume that the "objseq" is a sequence  
     #  (tuple or list) of objnotes. An "objnotes" consists of a sequence 
     #  pairing of an object and its annotation.  Example:
     #       objseq = [ (obj1, 'First thing'), (obj2, 'Second thing') ]
     #  Use an empty string like "" to explicitly blank out annotation.

     def inbatch( self, objseq, table=Base.tab0 ):
          '''Pickle and compress sequence of annotated objects; insert.'''
          self.createtable( table ) 
          #    ^ serves also to check table's existence.
          s  = "INSERT INTO %s " % table
          v  = "VALUES (null, strftime('%s','now'), ?, ?)"
          #                   ^SQLite's function for unix epoch time.
          sql = ' '.join([s, v])
          def generate_parlist():
               for i in objseq:
                    obj, notes = i
                    parlist  = [ notes, ysql.Binary(pzdumps(obj)) ]
                    yield parlist
                    #     ^ using generator for parameter list.
          self.proceed( sql, generate_parlist() ) 
          #    inserting 100,000 rows takes about 10 seconds.

     #  objseq can be generated on the fly. Just write a generator function, 
     #  and pass it along to pzgenerator [for illustration, see copy].

     def ingenerator( self, generate_objnotes, table=Base.tab0 ):
          '''Pickle and compress via generator function, then insert.'''
          #  generator should yield an objseq element like this: (obj, notes)
          self.inbatch( [x for x in generate_objnotes], table )

          #  TIP:  generate computationally intense results, then 
          #  pass them to ingenerator which will warehouse them. Instantly 
          #  access those pre-computed results later by subquery on notes.


     def insert( self, obj, notes='#0notes', table=Base.tab0 ):
          '''Pickle and compress single object; insert with annotation.'''
          self.inbatch( [(obj, notes)], table )

          #  CAVEAT: if you have *lots* of objects to insert individually 
          #  this repeatedly will be slow because it commits after every 
          #  insert which causes sqlite to sync the inserted data to disk. 
          #       REMEDY: prepare your annotated objects in objseq form, 
          #       then use inbatch or ingenerator instead.  



class Annex( Insertion ):
     '''_______________ Add macro-objects (files, URL content) to DATABASE'''

     def inweb( self, URL, notes='', table=Base.tab0 ):
          '''Pickle and compress URL content, then insert into table.'''
          #     put URL address in quotes including the http portion.
          if not notes:
               #  let notes be an empty string if you want the URL noted,
               #            else notes will be that supplied by argument. 
               notes = URL
          import urllib2
          webpage = urllib2.urlopen( URL )
          webobj  = webpage.read()
          self.insert( webobj, notes, table )

     def file2string( self, filename ):
          '''Convert any file, text or binary, into a string object.'''
          #          put filename in quotes including the path.
          f = open( filename, 'rb' )
          #                   'read binary' but also works here for text;
          #                    line-end conversions are suppressed.
          strobj = f.read()
          f.close()
          return strobj
          #      ^       Note: even if the file originally was iterable.
          #              that STRING OBJECT will NOT be iterable.
          #  If the file was text, one could regain iteration by writing
          #  back to a disk file, or by doing something in-memory, e.g. 
          #                      import cStringIO as yString
          #                      #   or  StringIO
          #                      inmemoryfile = yString.StringIO( strobj )
          #                      for line in inmemoryfile:
          #                          ... process the line ... 
          #                      inmemoryfile.close()
          #       Iteration is not necessary if the string object is 
          #       to be used in its entirety, e.g. some boilerplate text.
          #       Getting text snippets from a database is generally 
          #       faster than reading equivalent text files from disk.

     def infile( self, filename, notes='', table=Base.tab0 ):
          '''Pickle and compress any file, then insert contents into table.'''
          #    put filename in quotes including any needed path.
          #    file can be binary or text, including image.
          if not notes:
               #  let notes be an empty string if you want filename noted,
               #                 else notes will be that supplied by argument. 
               notes = filename
          self.insert( self.file2string( filename ), notes, table )


class Answer( Base ):
     '''_______________ Single item answer shouted out.'''
     #  see class Selection for reading and retrieval.

     def responder( self, klass, tupler, response ):
          '''(see def respond in abstract superclass Base.)'''
          #   different behavior depending on subclass...
          if klass == 'Answer':
               response[0] = tupler
               #        ^ we only expect a single answer.
          if klass == 'Subquery':
               kid, tunix, notes, pzblob  =  tupler
               obj = pzloads( pzblob )
               response[kid] = [ tunix, notes, obj ]
               #  each item in response DICTIONARY has a key kid          <= 
               #  (same as in the table), and it is a list consisting of  <= 
               #       timestamp, notes, and original object              <= 
               #                         (decompressed and unpickled).    <= 

     def shout( self, question, table=Base.tab0 ):
          '''Shout a question; get a short answer.'''
          sql = "SELECT ( %s ) FROM %s" % (question, table)
          response = self.respond( 'Answer', sql )
          return response[0][0]
          #      None, if the table is empty.

     def lastkid( self, table=Base.tab0 ):
          '''Get primary key ID of the last insert.'''
          maxkid = self.shout( "MAX( kid )", table ) 
          if maxkid == None:
               maxkid = 0
          return maxkid

     def lastsec( self, table=Base.tab0 ):
          '''Get time in unix seconds of the last insert.'''
          tmax = self.shout( "MAX( tunix )", table ) 
          if tmax == None:
               tmax = 0
          return tmax
          #      e.g.  1256856209
          #    seconds ^ elapsed since 00:00 UTC, 1 January 1970

          #       _____ Automatic timestamp ("tunix")
          #  As each object enters the database we also add an  
          #  epoch timestamp.  Rather than import another Python 
          #  module for that purpose, we use SQLite's (faster) 
          #  functions: strftime and datetime.

     def lastdate( self, table=Base.tab0 ):
          '''Get local date/time of the last insert.'''
          tmax = self.lastsec( table ) 
          q = "datetime( %s, 'unixepoch', 'localtime')" % tmax
          if not tmax:
               return ' :: lastdate not applicable.'
          else:
               return  self.shout( q, table )
               # e.g.  u'2009-10-29 15:43:29' 



class Util:
     '''_______________ Utility methods for keys, subquery, comma'''

     def reverse_dickeys( self, dictionary, recentfirst=True ):
          '''List dictionary keys: sorted reverse (chronological) order.'''
          dickeys = dictionary.keys()
          dickeys.sort()
          if recentfirst:
               dickeys.reverse()
          return dickeys


     #           _____ SUBQUERY regex style for LIKE  
     #
     #  Here we want the dictionary to consist of items which have
     #  notes containing " gold " (using LIKE):
     #
     #             dic = I.dicsub("WHERE notes LIKE '% gold %'")
     #
     #  Percent    % in LIKE is the regex equivalent of star      "*"
     #             ^wildcard for 0 or more characters
     #  Underscore _ in LIKE is the regex equivalent of period    "."
     #             ^single character
     #  Escape     ! in LIKE is the regex equivalent of backslash "\"
     #                  (or one can specify an escape character 
     #                   by appending, for example, "ESCAPE '\'"
     #                   at the end of the subquery.)
     #  [a-c]   same in LIKE for character ranges, exclude by     "^"

     #           _____ SUBQUERY regex style for GLOB 
     #
     #  For SQLite: the GLOB operator is similar to LIKE but uses the 
     #  Unix file globbing [dissimilar to grep] syntax for its wildcards. 
     #             ? in GLOB is the regex equivalent of period    "."
     # 
     #                  GLOB is case sensitive, unlike LIKE.            <= 
     #  But LIKE is case sensitive for unicode characters beyond ASCII.
     #  Both GLOB and LIKE may be preceded by the NOT keyword 
     #  to invert the sense of the test.

     def notesglob( self, parlist ):
          '''Create a CONJUNCTIVE subquery using GLOB with placeholder.'''
          #           ^i.e. each term in parlist is an "AND" search term.
          s = ""
          for i in parlist:
               s += "AND notes GLOB ? "
               #    use placeholder ^ rather than i for security!
          return s.replace('AND', 'WHERE', 1)
          #        replace only the first occurrence


     #  Everyone is going to have a unique style of writing out their notes 
     #  so that it can be efficiently searched. Tags are helpful, but they 
     #  are optional within notes. Tags are useful for creating indexes.
     #       Define a "TAG" to be text in notes prefixed by "#" 

     def comma2list( self, csvstr, wild=True ):
          '''Convert comma separated values within string to a parameter list 
               >>> #     !! white spaces within string are significant !! 
               >>> print comma2list('#paris, agent007 ,#scheme')
               ['*#paris*', '* agent007 *', '*#scheme*']

             Empty or single entry (without comma) csvstr is acceptable.
             Empty string '' will select everything if wild is True.
             Unlike official csv, internal quotes should not be used;
             for simplicity, comma itself is not meant to be escaped.
          '''
          #  wild will conveniently include stars on both ends...
          if wild:
               parlist = [ '*%s*' % i for i in csvstr.split(',') ]
          else:
               parlist = [  '%s'  % i for i in csvstr.split(',') ]
               #  manually add wildcards as needed to csvstr for faster regex.
          return parlist 
          #
          #  TIP: for faster execution, list most discriminating values first;
          #       use wild=False for exact search wherever possible.
          #
          #       To form csvstr out of string variables a, b, c:
          #               csvstr = ','.join( [a, b, c] )
          #       See function comma after Main class.



class Deletion( Base, Util ):
     '''_______________ Deletion methods; also used for queue POP'''

     def deletesub( self, subquery, parlist=[], table=Base.tab0 ):
          '''Delete row(s) matching the subquery.'''
          sql = 'DELETE FROM %s %s'  % ( table, subquery )
          #    use ? placeholder(s) for security^
          self.proceed( sql, [ parlist ] ) 

     def deletekid( self, kid, table=Base.tab0 ):
          '''Delete single row with primary key kid.'''
          subquery = 'WHERE kid = ?' 
          self.deletesub( subquery, [ kid ], table )

     def deletecomma( self, csvstr, table=Base.tab0, wild=True ):
          '''Delete row(s): notes match comma separated values in string.'''
          parlist = self.comma2list( csvstr, wild )
          self.deletesub( self.notesglob(parlist), parlist, table )

     def delete( self, dual, table=Base.tab0, wild=True ):
          '''Alias "delete":            deletekid OR deletecomma.'''
          #  assuming dual is either an ^integer  OR ^csvstr string...
          if isinstance( dual, int ):
               self.deletekid(   dual, table       )
          else:
               self.deletecomma( dual, table, wild )

     def droptable( self, table=Base.tab0 ):
          '''Delete a table: destroys its structure, indexes, data.'''
          sql = 'DROP TABLE %s' % table 
          try:
               self.proceed( sql ) 
          except:
               if DEBUG:
                    print " ?? droptable: no table to delete."
          return " :: droptable: done."

     #  Delete a SQLite database file like a normal file at OS level.



class Subquery( Util, Answer, Deletion ):
     '''_______________ SUBQUERY table, get dictionary. POP QUEUE.'''
     #  "kid" serves as key for both retrieved dictionaries and database.
     #
     #            _____ SUBQUERY : SECURITY NOTICE.
     #       Per the DB-API recommendations,
     #            subquery should use ? as parameter placeholder
     #            to prevent SQL injection attacks; see Endnotes.
     #       Such a placeholder does not take a clause, simply values.
     #            Obscured fact: table names cannot be parametized.
     #       parlist shall be the parameter list which sequentially
     #            corresponds to the placeholder(s).
     #            parlist should be empty [] if no placeholders are used.

     def dicsub(self, subquery='', parlist=[], table=Base.tab0, POP=False):
          '''Subquery table to get objects into response dictionary.'''
          a = 'SELECT kid, tunix, notes, pzblob FROM %s %s'
          sql = a % ( table, subquery )
          response = self.respond( 'Subquery', sql, parlist )
          if POP:
               self.deletesub( subquery, parlist, table )
          return response

     #       __________ Using POP for QUEUE purposes         ___ATTN___ 
     #
     #  After y_serial retrieves entities that match a subquery pattern, 
     #  it can optionally delete them.
     #  
     #       POP = True, herein means "Treat as queue"               <=!!! 
     #       whereby objects matching a subquery are 
     #       retrieved and then DELETED.
     #     
     #       POP = False, means "retrieve but DO NOT delete."        <=!

     def diclast( self, m=1, table=Base.tab0, POP=False ):
          '''Get dictionary with last m consecutive kids in table.'''
          kid = self.lastkid( table ) - m 
          return self.dicsub('WHERE kid > ?', [kid], table, POP )

     def diccomma( self, csvstr, table=Base.tab0, wild=True, POP=False ):
          '''Get dictionary where notes match comma separated values.'''
          parlist  = self.comma2list( csvstr, wild )
          subquery = self.notesglob( parlist )
          return self.dicsub( subquery, parlist, table, POP )

     def selectdic( self, dual=1, table=Base.tab0, POP=False ):
          '''Alias "selectdic":         diclast  OR diccomma.'''
          #  assuming dual is either an ^integer OR ^csvstr string...
          if isinstance( dual, int ):
               return self.diclast(  dual, table,       POP )
          else:
               wild = True
               #      ^constrained for dual usage
               return self.diccomma( dual, table, wild, POP )



class Display( Subquery ):
     '''_______________ View subquery via pretty print'''

     def viewsub(self, subquery='', parlist=[], table=Base.tab0, POP=False):
          '''Subquery, order keys, and print qualified dictionary.'''
          dic        = self.dicsub( subquery, parlist, table, POP )
          dickeylist = self.reverse_dickeys( dic )
          diclen     = len( dickeylist )
          print "\n ::  View in reverse chronological order :: "
          print   " :: ----------------------------------------"
          #        ^^^^ "grep -v '^ :: '"  to get rid of labels.
          if diclen:
               for kid in dickeylist:
                    [ tunix, notes, obj ] = dic[ kid ]
                    print " ::   kid:  %s   (%s secs)" % (kid, tunix)
                    print " :: notes: ", notes
                    if isinstance( obj, str ):
                         if len( obj ) > 1000:
                              end = '[... Display limited ...]'
                              obj = ' '.join( [obj[0:1000], end] )
                         print obj
                    else:
                         try:
                              pprint.pprint( obj )
                         except:
                              print " !! Display: object not printable."
                    print " :: ----------------------------------------"
               print " :: Display: MATCHED  %s objects." % diclen
               if POP:
                    print "             and POP deleted them." 
          else:
               print " !!  Display: NOTHING matched subquery !!"

     def viewlast( self, m=1, table=Base.tab0, POP=False ):
          '''Print last m consecutive kids in table.'''
          kid = self.lastkid( table ) - m 
          self.viewsub('WHERE kid > ?', [kid], table, POP )

     def viewcomma(self, csvstr='', table=Base.tab0, wild=True, POP=False):
          '''Print dictionary where notes match comma separated values.'''
          parlist  = self.comma2list( csvstr, wild )
          subquery = self.notesglob( parlist )
          self.viewsub( subquery, parlist, table, POP )

     def view( self, dual=1, table=Base.tab0, POP=False ):
          '''Alias "view":              viewlast OR viewcomma.'''
          #  assuming dual is either an ^integer OR ^csvstr string...
          if isinstance( dual, int ):
               self.viewlast(  dual, table,       POP )
          else:
               wild = True
               #      ^constrained for dual usage
               self.viewcomma( dual, table, wild, POP )
          return ''



class Latest( Display ):
     '''_______________ Retrieve the latest qualified object "omax" '''
     #  TIP :: snippets from omaxsub could be helpful in one's program.

     def omaxsub(self, subquery='', parlist=[], table=Base.tab0, POP=False):
          '''Get the latest object omax which matches subquery.'''
          dic = self.dicsub( subquery, parlist, table )
          dickeylist = self.reverse_dickeys( dic )
          diclen  = len( dickeylist )
          #         count how many matched subquery
          if diclen:
               keymax = dickeylist[0]
               # ^MAX DIC KEY of the LATEST matching subquery.
               #        dic[keymax][0] corresponds to tunix.
               #        dic[keymax][1] corresponds to notes.
               omax   = dic[keymax][2]
               #  ^this is the LATEST OBJECT matching subquery.
               if POP:
               #    ^queue-like deletion of only single object:
                    self.deletekid( keymax, table )
          else:
               omax = None
          return omax

     def omaxlast( self, n=0, table=Base.tab0, POP=False ):
          '''Most quickly get the latest n-th object using key index.'''
          #               n = 0,1,2,...  assuming consecutive kids.
          #  Avoiding LIKE or GLOB enhances performance in large tables.
          #  Also we are not reading an entire table into memory since a
          #  blank subquery would have put an entire table into dictionary.
          #
          subquery = "WHERE kid=(SELECT MAX(kid) - ? FROM %s)" % table
          obj = self.omaxsub( subquery, [n], table, POP )
          if DEBUG and obj == None:
               print " !! omaxlast: that kid does not exist."
          return obj

     def omaxcomma( self, csvstr, table=Base.tab0, wild=True, POP=False ):
          '''Get latest object where notes match comma separated values.'''
          parlist  = self.comma2list( csvstr, wild )
          subquery = self.notesglob( parlist )
          return self.omaxsub( subquery, parlist, table, POP )

     def select( self, dual=0, table=Base.tab0, POP=False ):
          '''Alias "select":            omaxlast OR omaxcomma.'''
          #  assuming dual is either an ^integer OR ^csvstr string...
          if isinstance( dual, int ):
               obj  = self.omaxlast(  dual, table,       POP )
          else:
               wild = True
               #      ^constrained for dual usage
               obj  = self.omaxcomma( dual, table, wild, POP )
          return obj

     def getkid( self, kid, table=Base.tab0, POP=False ):
          '''Retrieve a row given primary key kid, POP optional.'''
          subquery = 'WHERE kid = ?' 
          return self.omaxsub( subquery, [ kid ], table, POP )
          #  added 2010-06-07  can't believe this was missing ;-)



class Oldest( Latest ):
     '''_______________ Retrieve the oldest qualified object "omin" '''

     def ominfirst( self, n=0, table=Base.tab0, POP=False ):
          '''Most quickly get the oldest n-th object using key index.'''
          #               n = 0,1,2,...  assuming consecutive kids.
          subquery = "WHERE kid=(SELECT MIN(kid) + ? FROM %s)" % table
          obj = self.omaxsub( subquery, [n], table, POP )
          #          ^works because its dictionary is single entry.
          if DEBUG and obj == None:
               print " !! ominfirst: that kid does not exist."
          return obj

     def fifo( self, table=Base.tab0 ):
          '''FIFO queue: return oldest object, then POP (delete) it.'''
          n   = 0
          POP = True
          return self.ominfirst( n, table, POP )



class Care( Answer, Deletion ):
     '''_______________ Maintenance methods'''

     def freshen( self, freshdays=None, table=Base.tab0 ):
          '''Delete rows in table over freshdays-old since last insert.'''
          #         freshdays could be fractional days, e.g. 2.501 days;
          #                   if it is None, then it's infinity.
          #                   if it is 0, then nothing will remain.
          if freshdays != None :
               max_tunix  = self.lastsec( table )
               freshsecs  = int( freshdays * 86400 )
               expiration = max_tunix - freshsecs
               sql        = "WHERE tunix <= ?"
               self.deletesub( sql, [expiration], table )

     def vacuum( self ):
          '''Defrag entire database, i.e. all tables therein.'''
          self.proceed( 'VACUUM' ) 

          #    _____ why VACUUM?
          #  "When an entity (table, index, or trigger) is dropped from 
          #  the database, it leaves behind empty space. This empty space 
          #  will be reused the next time new information is added; but 
          #  the database file might be larger than strictly necessary. 
          #  Also, frequent inserts, updates, and deletes can cause the 
          #  information in the database to become fragmented. The VACUUM 
          #  command cleans the main database by copying its contents
          #  to a temporary database file and reloading the original 
          #  database file from the copy. This eliminates free pages, 
          #  aligns table data to be contiguous, and otherwise cleans up 
          #  the database file structure." -- sqlite.org
          #  N.B. -  Surprising how much file size will shrink.

     def clean( self, freshdays=None, table=Base.tab0 ):
          '''Delete stale rows after freshdays; vacuum/defrag database.'''
          self.freshen( freshdays, table )
          self.vacuum()
          return ''



class Main( Annex, Oldest, Care ):
     '''_______________ Summary for use of a single database.'''
     pass
     #                  Base
     #        Insertion(Base)
     #  Annex(Insertion)
     #                                                              Util
     #                                       Answer(          Base)
     #                                               Deletion(Base, Util)
     #                        Subquery(Util, Answer, Deletion)
     #                Display(Subquery)
     #         Latest(Display)
     #  Oldest(Latest)
     #                                  Care(Answer, Deletion)



#  _______________ COPY functions (demonstration outside of Main class)
#                       also note how ingenerator is employed usefully.

def copysub( subquery, parlist, tablex, tabley, dbx=Base.db0, dby=Base.db0 ):
     '''Subselect from tablex, then copy to tabley (in another database).'''
     #          assume tablex is in dbx, and tabley is in dby;
     #          copying from *x to *y.
     if (tablex != tabley) or (dbx != dby):
          X = Main( dbx )
          dic = X.dicsub( subquery, parlist, tablex )
          dickeylist = X.reverse_dickeys( dic, recentfirst=False )
          #                     order keys chronologically ^
          #                     to preserve inserted ordering on copy.
          diclen  = len( dickeylist )
          #         count how many matched subquery
          if diclen:
               Y = Main( dby )
               def generate_objnotes():
                    for i in dickeylist:
                         notes = dic[i][1]
                         obj   = dic[i][2]
                         yield (obj, notes)
               Y.ingenerator( generate_objnotes(), tabley )
               #   generator copies objects & notes from qualified dictionary.
               #   Timestamps are fresh, i.e. not preserved from old table.
               if DEBUG:
                    p = ( diclen, tablex, tabley )
                    print " :: copysub:  %s objects from %s to %s." % p
          else:
               if DEBUG:
                    p = (         tablex, tabley )
                    print " !! copysub:  NOTHING from %s to %s." % p
     else:
          print " !! copysub: table or database name(s) must differ."

def copylast( m, tablex, tabley, dbx=Base.db0, dby=Base.db0 ):
     '''Copy last m consecutive kids in tablex over to tabley.'''
     A = Answer()
     kid = A.lastkid( tablex ) - m 
     copysub('WHERE kid > ?', [kid], tablex, tabley, dbx, dby )

def comma( *string_args ):
     '''Join string-type arguments with comma (cf. csvstr, comma2list).'''
     #                   ^which may include regular expression...  Essential <=!
     return ','.join( string_args )

def copycomma( csvstr, tablex, tabley, dbx=Base.db0, dby=Base.db0, wild=True ):
     '''Subselect by comma separated values from tablex, then copy to tabley.'''
     U = Util()
     parlist  = U.comma2list( csvstr, wild )
     subquery = U.notesglob( parlist )
     copysub( subquery, parlist, tablex, tabley, dbx, dby )

def copy( dual, tablex, tabley, dbx=Base.db0, dby=Base.db0, wild=True ):
     '''Alias "copy":              copylast OR copycomma'''
     #  assuming dual is either an ^integer OR ^csvstr string...
     if isinstance( dual, int ):
          copylast(  dual, tablex, tabley, dbx, dby       )
     else:
          copycomma( dual, tablex, tabley, dbx, dby, wild )



# ============================= BETA ================================================= 

import random

class Farm:
     '''_______________ Start a farm of databases for concurrency and scale.'''
     #  (Dependencies: Insertion, Deletion classes; copy function.)
     #  Reduces probability of locked database, prior to insert, by 
     #  writing to a random "barn" database within "farm" directory.  
     #  A barn is intended for temporary storage until harvest at  
     #  which time the objects are moved to some central database. 
     # 
     #  The number of barns increases concurrent write possibilities.
     #  The variable "size" determines frequency of harvest.
     #  The two variables should be optimized for your situation.

     dir0   = '/var/tmp/yaya/db/y_serial_farm'
     barns0 = 9
     #        ^ can be greatly increased without any tax on memory,
     #          for it only increases the number of files in dir.

     #  To set maxbarns, estimate the number of writes per second, 
     #  and divide that by, say, 4.

     def __init__( self, dir=dir0, maxbarns=barns0 ):
          '''Set directory for the farm of maxbarns database files.'''
          #  Be sure it exists at OS level with appropriate permissions.
          self.dir = dir
          if not self.dir.endswith( '/'):
               self.dir += '/'
          self.maxbarns = maxbarns

     def barn( self, n ):
          '''Prepend directory to numbered database filename.'''
          return ''.join( [self.dir, "barn%s.sqlite" % n ] )
          #  e.g.   /var/tmp/yaya/db/y_serial_farm/barn7.sqlite

     def farmin( self, obj, notes='#0notes', table=Base.tab0, n=0 ):
          '''Insert an object with notes into barn(n).'''
          #  n should be random in range(maxbarns), see plant method.
          I = Insertion( self.barn(n) )
          I.insert( obj, notes, table )

     def reap( self, dual, tablex, tabley, n, dby=Base.db0, wild=True ):
          '''Move dual under tablex in barn(n) to tabley in dby.'''
          try:
               copy( dual, tablex, tabley, self.barn(n), dby, wild )
               D = Deletion( self.barn(n) )
               D.delete( dual, tablex, wild )
          except:
               if DEBUG:
                    print " :: reap: skipped barn %s " % n

     #  Objects enter a barn singularly via farmin, however, they move out 
     #  rapidly by generator via reap.  The harvest method spells out the 
     #  stochastic condition for reap to occur (after each farmin insert). 
     #  This technique avoids frequent and expensive query of barn contents.
     #  * See plant method below for typical usage.

     def harvest(self, dual, tablex,tabley, n, dby=Base.db0, wild=True, size=10):
          '''After farmin, reap dual under tablex barn(n) by size expectation.'''
          #  Prioritize:  decreasing size increases the probability of reap.
          #  If dual='' and we harvest after every farmin, we can statistically  
          #  expect size number of objects moved whenever reap is triggered.
          #  (For size=10, reap may move from 1 to 47 objects, 10 on average.)
          if size * random.random() < 1:
               self.reap( dual, tablex, tabley, n, dby, wild )
          else:
               if DEBUG:
                    print " :: harvest: nothing from barn%s" % n

     #  TIP: use reap to flush remaining objects at the close your script.

     def cleanfarm( self, freshdays=None, table=Base.tab0 ):
          '''Delete stale rows after freshdays; vacuum/defrag barns.'''
          for n in range( self.maxbarns ):
               try:
                    C = Care( self.barn(n) )
                    C.clean( freshdays, table )
               except:
                    if DEBUG:
                         print " :: cleanfarm: skipped barn %s" % n
          if DEBUG:
               print " :: cleanfarm: VACUUMed barns in %s" % self.dir

     def plant( self, obj, notes='#0notes', table=Base.tab0, dby=Base.db0 ):
          '''FARM SUMMARY: farmin insert with generic self-cleaning harvest.'''
          size = 10
          wild = True
          #      ^constrains, also: dual='' and only one table name.
          if obj == 'reap_ALL_BARNS':
          #  A bit odd, but it beats writing out the iteration later...    ;-)
               notes = ' :: plant:  object and notes were not inserted.'
               for n in range( self.maxbarns ):
                    self.reap( '', table, table, n, dby, wild )
          else:
               n = random.randrange( self.maxbarns )
               #   why random? to minimize conflicts other farmin operations.
               self.farmin( obj, notes, table, n )
               #    ^inserts obj into table in some random barn(n).
               #  print "farmin barn%s" % n
               self.harvest( '', table, table, n, dby, wild, size )
               #   ^harvest whenever around size accumulates in a random barn.
          if 100000 * random.random() < 1:
          #  vacuum of all barns approximately every 100,000 inserts.
               self.cleanfarm()



# =========================== ENDNOTES =============================================== 

#  What problem does y_serial solve beyond serialization?
#  
#       "pickle reads and writes file objects, it does not handle the 
#       issue of naming persistent objects, nor the (even more complicated) 
#       issue of concurrent access to persistent objects. The pickle module 
#       can transform a complex object into a byte stream and it can 
#       transform the byte stream into an object with the same internal 
#       structure. Perhaps the most obvious thing to do with these byte
#       streams is to write them onto a file, but it is also conceivable 
#       to send them across a network or store them in a database." 
#            http://docs.python.org/library/pickle.html


#  y_serial takes a couple of minutes to write a MILLION annotated objects  
#           (that includes serialization and compression),
#           which will consume at least 37MB for a single sqlite file, 
#           on a rather antiquated 32-bit commodity desktop, 
#           blink of an eye for read access using GLOB regex.
#           Generally much faster than comparable DB-API use with PostgreSQL.


#  "To BLOB or Not To BLOB: Large Object Storage in a Database or a Filesystem?
#                           by Russell Sears; Catharine Van Ingen; Jim Gray
#       Paper submitted on 26 Jan 2007 to http://arxiv.org/abs/cs.DB/0701168 
#  
#  This paper looks at the question of fragmentation [...] objects smaller than
#  256KB are best stored in a database while objects larger than 1M are best 
#  stored in the filesystem. Between 256KB and 1MB, the read:write ratio and 
#  rate of object overwrite or replacement are important factors." 
#  
#       Generally speaking, database queries are faster than file opens,  
#       however, filesystems are optimized for streaming large objects.
#       That paper shows how important it is to keep a lean database, 
#       which is why we wrote the sqliteclean function. 
#
#       Objects under 1MB will do fine since they will be compressed 
#       before insertion into the database. If you have many larger objects, 
#       we included the alternative gzip compressed file solution.


#   SUPPLEMENT
# _______________ pz FUNCTIONS for FILE.gz        [not database related]
#                 Using compression with pickling.
#                 Source: recipe 7.3,  Python Cookbook, second edition.
#
     # ___ATTN___ Individual pickled items are already compressed by above.
     #            Below we place those items in a file which is gzipped.
     #            So there are two compression stages...
     #            the second may not squeeze out much, but we might as
     #            well gzip as long as we are writing to a file.

import gzip
#      ^compression for files.

def pzdump(filename, *objects):
     '''Pickle and zlib-compress objects, then save them in a gz file.'''
     fil = gzip.open(filename, 'wb')
     for obj in objects: 
          yPickle.dump( pzdumps(obj), fil, pickle_protocol)
     fil.close()

     #  The protocol is recorded in the file together with the data, so 
     #  Pickler.load can figure it out. Just pass it an instance of a file 
     #  or pseudo-file object with a read method, and Pickler.load returns 
     #  each object that was pickled to the file, one after the other, 
     #  and raises EOFError when the file's done. We wrap a generator 
     #  around Pickler.load, so you can simply loop over all recovered 
     #  objects with a for statement, or, depending on what you need, 
     #  you can use some call such as list(pzload('somefile.gz')) 
     #  to get a list with all recovered objects as its items.

def pzload(filename):
     '''Iterate zlib-compressed pickled objects from a gz file.'''
     fil = gzip.open(filename, 'rb')
     while True:
          try: 
               yield   yPickle.load(fil) 
          #          ^hang on to the compressed version for now,
          #                   decompress later as needed.
          #    ^ iterator
          except EOFError: 
               break
     fil.close()

     #  Example of iteration use...
     #    for i in pzload(filename):  print pzloads(i)
     #    #   each item gets printed after decompression.

def pzlist(filename):
     '''List of zlib-compressed pickled objects from a gz file.'''
     return list( pzload( filename ) )

def oblist(filename):
     '''List of zlib-decompressed pickled objects from a gz file.'''
     return [ pzloads(i) for i in pzload(filename) ]

     #  Another example of iteration use, ASSIGN VARIABLES...
     #    [x, y, z] = oblist(filename)
     #    #           assuming three items in the file.


#  ================== DATABASE versus FILE.gz ======================= 
#
#       Putting all pz objects into a file would be suitable 
#       where the collection of such objects is fairly static 
#       and not so large in quantity. Use a database like 
#       SQLite or PostgreSQL if:
#            - the situation is dynamic, i.e. pz objects need
#                 to be appended or deleted often.
#            - particular pz objects are needed
#                 (using files, all pz objects have to be
#                 unpacked, then picked over).
#            - the database can annotate or index the contents of 
#                 of pz objects (use SQL to then cherry-pick).
#            - the objects are generally under 1MB; see Endnotes.
#
#  # _______________ 2009-08-24  warehouse objects in a file.gz
#       import y_serial
#       fname = '/tmp/y_serial.gz'
#  
#       item = 'This is a string to TEST some file.gz'
#       y_serial.pzdump( fname, item, item, item, item, item )
#       print y_serial.oblist( fname )
#
#  ================================================================== 



#  Why BLOB dictionaries for storing schema-less data?
#  
#       * Dictionaries are very suitable objects 
#         which are variable length, arbitrarily nestable, 
#         and can contain arbitrary objects...
#
#       Bret Taylor wrote a wonderful post about how Friendfeed uses 
#       MySQL to store schema-less data, 
#       http://bret.appspot.com/entry/how-friendfeed-uses-mysql 
#       which got me thinking about the details and its use for 
#       any Python program.  (Thanks very much, Bret! ;-)
#
#       Friendfeed uses a database with schema-less data, where 
#       dictionaries are zlib-pickled and then inserted into MySQL. 
#       Other tables then index the primary. Easier to shard, and 
#       avoids JOINs -- thus conceptually, python code on the 
#       dictionary structure replaces SQL code. The database 
#       merely becomes a fancy hash table with fast access.
#       Nice for rapid development because dictionaries are 
#       easily modifiable and themselves have fast key access.
#       Better for database load distribution and maintenance, 
#       plus it avoids scary table conversions requiring downtime.


#  "Stop calling me NoSQL" by Dhananjay Nene
#            http://blog.dhananjaynene.com/2009/10/stop-calling-me-nosql/
# 
#       You see unlike RDBMS, I don't require that data be clearly 
#       split into tables, columns and rows. I can work with data the way 
#       it is most naturally represented:  as a tree of individual data 
#       fields, lists, arrays, dictionaries, etc. Also I do not require 
#       that you always clearly define each and every possible schema
#       element before being able to store data corresponding to the 
#       schema. I can happily accept a schema dynamically or even 
#       work without a schema. Some of my early forms were based on 
#       key value pairs stored as B-Trees (eg. Berkeley DB).  Over the 
#       years people have figured out ways to represent the data as 
#       a set of decomposed document elements, store data spread across 
#       a cluster, replicate it for better availability and fault tolerance, 
#       and even perform post storage processing tasks using map-reduce 
#       sequences.  But really what separates me from my cousin and other 
#       storage systems is that I don't make demands on the data -- I take 
#       it in its naturally found form and then store it, replicate it, 
#       slice it, dice it and glean information out of it. And therein 
#       lies my true identity -- I will work with data the way the data 
#       is best represented with all its arbitrary inconsistencies and 
#       inabilities to always clearly specify a constraining schema.



#  What can the PICKLE/cPickle module store? and what about json?
#  
#      * All the native datatypes that Python supports: booleans, integers, 
#           floating point numbers, complex numbers, strings, bytes objects, 
#           byte arrays, and None.
#      * Lists, tuples, dictionaries, and sets containing 
#           any combination of native datatypes.
#      * Lists, tuples, dictionaries, and sets containing any combination of 
#           lists, tuples, dictionaries, and sets containing any combination 
#           of native datatypes (and so on, to the maximum nesting level 
#           that Python supports).
#      * Functions, classes, and instances of classes (with CAVEATS):
#           pickle can save and restore class instances transparently, 
#           however the class definition must be importable and live in the 
#           same module as when the object was stored.  picklable functions 
#           and classes must be defined in the top level of a module.
#                [ Most likely reason why pzget gets CHOKED. 
#                  (Hack: insert the defining text, then exec it later.) ]
#  
#      Good reference: http://diveintopython3.org/serializing.html
#      Also includes a comparative review of the json module 
#      introduced as of Python v2.6 -- which is text-based serialization.
#
#           _____ json versus pickle 
#            Few reasons why we opted for pickle instead of json:
#            - human-readability is not a primary concern
#                 since the database could care less.
#            - json does not distinguish between tuples and lists.
#            - json cannot handle complex Python objects 
#                 without additional en/decoding.
#                 (and why worry about internal structures?)
#            - since json uses utf-8, this may fail in some cases:
#                 obj == json.loads(json.dumps(obj))
#            - we are not handing off the serialized item 
#                 to be read by another language.
#            - as for SECURITY, we are not accepting any serialized 
#                 item from an untrusted source into the database.
#                 y_serial's particular use of pickle is safe.

#            _____ short digression on pickle security risk
#
#  Generally, never unpickle an untrusted string whose origins are dubious, 
#  e.g. strings read from a socket or public webpage. So should one 
#  sanitize and encrypt such strings before the pickle stage? 
#  No, that would not be necessary.
#  
#       pickle uses a simple stack language that allows the creation 
#       of arbitrary python structures, and execute them. This stack 
#       language allows you to import modules (the 'c' symbol), and 
#       apply arguments to callables (the 'R' symbol), thus causing code 
#       to be run. Combine this with the python built-in methods eval 
#       and compile and you have the perfect vehicle for an 
#       "unpickle attack." 
#
#       For more details, also see excellent article by Nadia Alramli, 
#       http://nadiana.com/python-pickle-insecure
#            [And thanks so much, Nadia, for personally clarifying 
#             the difference between untrusted text string and 
#             the string derived from pickling such a thing.]
#
#  Some naive methods proposed to "encrypt" strings before pickle:
#  
#       import string
#       rot13key = string.maketrans(
#               'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz', 
#               'NOPQRSTUVWXYZABCDEFGHIJKLMnopqrstuvwxyzabcdefghijklm')
#       def rot13( text ):
#            '''Because Python v3 will discontinue .encode('rot13')'''
#            return text.translate( rot13key )
#            #  notice that, text == rot13(rot13( text ))
#  
#       #  Also look at: .encode('base64') with .decode('base64') 
#
#  But such effort is misplaced, because the security risk is not pickling 
#  untrusted data, but rather unpickling untrusted pickles.  The malicious
#  person must be in a position to be able to modify the pickle or replace 
#  it somehow (e.g. over a communication channel that is not secure). 
#  Pickles over public channels should be crytographically signed;
#  N. Alramli cites http://mumrah.net/2009/09/making-pythons-pickle-safer/
#
#  Never load untrusted pickle streams. Obey this frequently cited warning: 
#  "Never unpickle data received from an untrusted or unauthenticated 
#  source." 
#              y_serial complies, so don't worry, it's safe ...    
#                       (just keep your database secure from tampering ;-)


#  What does Google use in-house? not json, but rather "PROTOCOL BUFFERS," 
#  http://code.google.com/apis/protocolbuffers/docs/pythontutorial.html
#       Very nice for multi-lingual serialization, e.g. with C++ 
#       and Java, but its focus is schema-less messages, 
#       not Python objects.


#  What about YAML?
#       Saw http://pyyaml.org/wiki/PyYAML, which apparently 
#       offers "high-level API for serializing and deserializing 
#       native Python objects." In short, we did not have the time 
#       to test its reliability. Speed is relatively slow.


#  In summary, SERIALIZATION methods which are human readable, 
#  e.g. json or YAML, are much more slower than cPickle. 
#  If readability by other languages, e.g. C++ or Java, is 
#  not a necessary requirement, cPickle takes the cake. 


#  SQLite "maximum" number of TABLES
#
#  "The more tables you have, the slower the first query will run
#  and the more memory SQLite will use.  For long-running applications
#  where the startup time is not a significant factor, 100s or
#  1000s of tables is fine.  For a CGI script that starts itself
#  up anew several times per second, then you should try to keep
#  the number of tables below a 100, I think.  -- D. Richard Hipp"
#     http://www.mail-archive.com/sqlite-users@sqlite.org/msg14057.html


#  SQLite INTEGER PRIMARY KEY (regarding our kid)
#
#  For performance reasons we did not use the AUTOINCREMENT option 
#  which guarantees kid uniqueness, thus: "the largest ROWID is equal to 
#  the largest possible integer 9223372036854775807 then the database 
#  engine starts picking candidate ROWIDs at random until it finds 
#  one that is not previously used." -- That's a lot of rows!
#
#  Interesting to note that under the hood, "kid" as INTEGER PRIMARY KEY 
#  is just an alias for special column names ROWID, _ROWID_, or OID.


#  SQLite will work great as the database engine for low to medium traffic
#  WEBSITES (which is to say, 99.9% of all websites). Any site that gets 
#  fewer than 100K hits/day should work fine with SQLite.  SQLite has been 
#  demonstrated to work with 10 times that amount of traffic.
#  
#       see SQLite Appropriate Uses : http://www.sqlite.org/whentouse.html


#  SQLite SUBQUERY limits
#
#       The length of the LIKE or GLOB pattern is limited within SQLite to 
#       SQLITE_MAX_LIKE_PATTERN_LENGTH bytes. The default value of this 
#       limit is 50000. A typical computer can evaluate even a pathological 
#       LIKE or GLOB pattern of such size relatively quickly. 
#
#       Tip: the "notes" field could contain TAGS using the hashmark #.
#            That would focus searching for keywords.


#  Python Cookbook, second ed. #7.12. BLOB in SQLite does NOT work:
#
#       "The PySQLite Python extension offers function sqlite.encode 
#       to let you insert binary strings in SQLite databases." 
#       That function has been superceded by sqlite.Binary
#       And now ? replaces %s for more secure syntax.


#  DB-API string format v. SQL injection attacks
#  
#  2009-08-29 comp.lang.python.  Tim Cross illustrates why it's so important 
#  to use the DB API's own escaping functions. 
#  
#  >>  "somestring %s" % "foo" will work.
#  
#  > BAD IDEA when assembling SQL, unless you _like_ SQL-injection attacks:
#  > 
#  >      sql = "select * from users where name='%s' and password='%s'"
#  > 
#  >      # get some values from an untrusted user:
#  >      name = "administrator"
#  >      password = "' or 1=1; drop table users; --"
#  > 
#  >      cursor.execute(sql % (name, password))
#  >      # uh-oh!
#  
#  Of course, that OR 1=1 attack is at the publicly exposed prompt...
#  that clause attaches to the WHERE subquery which is always evaluates true!
#  oh, crap -- never knew how easy it was.
#  http://groups.google.com/group/comp.lang.python/browse_thread/thread/5fdaf7d1b46e6699


#      __________ PARAMETER SUBSTITUTION for sqlite3
#  
#  First how does one find out the substitution style per the DB-API?
#       >>> import sqlite3
#       >>> print sqlite3.paramstyle
#       qmark
#  
#  That means that a question mark "?" is used for SQL parameter substitions, and 
#  the second argument to the execute method is a *sequence* i.e. either a list 
#  or a tuple.  In the latter case, a typical gotcha error is using 
#       (myvariable) instead of (myvariable,) 
#  for a single variable substitution. Thus using a list is easier: [myvariable]
#
#  paramstyle can vary, for example, it's 'pyformat' for psycopg (PostgreSQL).
#  The string attribute paramstyle is apparently read-only.


#           _____ GOTCHA!  Parameter substitution for sqlite3
#  
#  After hours of insanity, I find out this very, very obscure fact:
#                           table names *cannot be parametized* !!
#  Indeed, one could easily think that a table name was a parameter, 
#  but then an injection attempt should not have access to other tables.
#  So be sure to never expose the table variable in public apps.
#
#  It also appears that under strict paramstyle, a placeholder cannot 
#  substitute a WHERE clause. This is very understandable since that's 
#  what a SQL injection attack wants to exploit. 


#  Alex Martelli: "I find the named-parameter binding style much more readable 
#  -- and sqlite3 supports it:
#  
#       c.execute('SELECT * FROM atable WHERE newID=:t', locals())
#  
#  Note: passing {'t': t} or dict(t=t) instead of locals() would be more 
#  punctiliously correct, but in my opinion it would interfere with readability 
#  when there are several parameters and/or longer names."


#      __________ Batch processing (re inbatch)
#
#  SQLite does fsync() 3 times per transaction to guarantee data integrity. 
#  So batch statements update the database in transactions 
#  (BEGIN TRANSACTION; ... COMMIT TRANSACTION;). Only 3 fsync are 
#  required per transaction, not per statement, and one also gets 
#  multi-statement atomicity, so all the changes make it to disk or 
#  none does.  cf. Firefox3 performance hit:
#       http://shaver.off.net/diary/2008/05/25/fsyncers-and-curveballs/
#
#  Richard Hipp concurs: "Using a transaction is the fastest way to 
#  update data in SQLite. After each transaction the SQLite engine 
#  closes and opens the database file. When SQLite opens a database file 
#  it populates the SQLite internal structures, which takes time. 
#  So if you have 100 updates and don't use a transaction then SQlite 
#  will open and close the database 100 times. Using transactions 
#  improves speed. Use them."
#
#  Ok, that said, it would seem to help if we did the following:
#            cur.execute( "BEGIN TRANSACTION" )
#               ... update stuff ...
#            cur.execute( "COMMIT TRANSACTION" )
#  but one should NOT do so, because the sqlite3 module has implicitly 
#  already taken care of this issue when it connects and con.commit()
#  
#       By default, sqlite3 opens transactions implicitly 
#       before a DML statement (INSERT/UPDATE/DELETE/REPLACE), 
#       and commits transactions implicitly before anything other 
#       than (SELECT/INSERT/UPDATE/DELETE/REPLACE).
#  
#       You can control the kind of "BEGIN" statements via the 
#       isolation_level parameter to the connect call,
#       or via the isolation_level property of connections.
#  
#       If you want autocommit mode, then set isolation_level 
#       to None [which does NOT begin transaction].  Otherwise 
#       the default will result in a plain "BEGIN" statement. 
#       One could also set it to an isolation level 
#       supported by SQLite: DEFERRED, IMMEDIATE or EXCLUSIVE.
#
#  y_serial uses IMMEDIATE; the differences are explained here:  
#       http://www.sqlite.org/lang_transaction.html 


#            _____ ENCODE TEXT in UTF-8
#  
#  Gerhard, back in 2007 said about the precursor to sqlite3: 
#       "SQLite databases store text in UTF-8 encoding. If you use pysqlite, 
#       and always use unicode strings, you will never have any problems with 
#       that. pysqlite does not rap on your knuckles if you store arbitrary 
#       encodings in the database, but you will feel sorry once you try to 
#       fetch the data."
#  For y_serial, this does not pertain to the object themselves because 
#  they are BLOBs, but it's relevant to the attached annotation notes.
#
#  [sqlite3 module will return Unicode objects for TEXT. If you wanted 
#   to return bytestrings instead, you could set con.text_factory to str.]


#            _____ converting timestamp ("tunix") into human form
#  
#  We rely on SQLite for time functions but in your own program you may 
#  find these python functions useful to convert unix epoch seconds:
#  
#            import time
#            #            e.g. let tunix = 1254458728
#  
#            def tunixdate( tunix ):
#                 return time.localtime(ticks)[:3]
#                 # e.g. (2009, 10, 1)
#  
#            def tunixclock( tunix ):
#                 return time.localtime(ticks)[3:6]
#                 # e.g. (21, 45, 28)
#  
#            def tunixtuple( tunix ):
#                 return time.localtime(ticks)[:6]
#                 # e.g. (2009, 10, 1, 21, 45, 28)


# ================================== TESTER ========================================== 

def tester( database=Base.db0 ):
     '''Test class Main for bugs. Include path for database file.'''
     ipass = 0
     print "Starting tester()...    for debugging y_serial module."
     if not DEBUG:
          print "[DEBUG switch is currently False.]"
     print "[Note: specify default database via db0 in class Base.]"
     print "Creating instance using database..."
     I = Main( database )
     print "  using database:", database
     ipass += 1
     I.createtable( 'ytest' )
     print "   created table: ytest"
     ipass += 1
     #   I.droptable( 'ytest' )
     # ^ comment out line to test a brand new table.
     # ================================================================== 
     print "INSERTING:  5 objects..."
     def generate_testitems( n ):
          for i in range( n ):
               objnotes = ( i, "testitem-%s" % i )
               yield objnotes
               #   ^yield, not return, for generators.
     I.ingenerator( generate_testitems(2), 'ytest' )
     ipass += 1
     #    --------------------------------
     tmp1 = ("Part of 3-tuple.", 98, 'Encode text in UTF-8.' )
     tmp2 = { 'spam' : 2 , 'eggs' : 43 }
     tmp3 = 'I aspire to be stringy.'
     I.inbatch([(tmp1, 'test #tuple'), (tmp2, 'test dictionary')], 'ytest')
     ipass += 1
     I.insert(tmp3, 'random string', 'ytest')
     ipass += 1
     #    --------------------------------
     lsec = I.lastsec( 'ytest' )
     print "     Checking epoch second of last insert: ", lsec
     ipass += 1
     ldate = I.lastdate( 'ytest' )
     print "     Checking local date/time of last insert: ", ldate
     ipass += 1
     lkid = I.lastkid( 'ytest' )
     print "     Checking last kid PRIMARY KEY: ", lkid
     ipass += 1
     #  print "     (Note: delete* methods v0.50 have passed inspection.)"
     # ================================================================== 
     print "     (Inserted and selected objects should be equivalent.)"
     print "     Trying omaxsub ..."
     got2 = I.omaxsub("WHERE notes GLOB ?", ['random*'], 'ytest')
     if got2 == tmp3:
          print "passed test: subquery."
          ipass += 1
     else:
          print "TEST FAIL!   subquery."
     #    --------------------------------
     print "     Trying omaxlast via .select ..."
     got1 = I.select( 1, 'ytest' )
     if got1['eggs'] == tmp2['eggs']:
          print "passed test: seek key ID."
          ipass += 1
     else:
          print "TEST FAIL!   seek key ID."
     #    --------------------------------
     print "     Trying omaxcomma ..."
     got3 = I.omaxcomma( comma('rand*','*ring'), 'ytest' )
     if got3 == tmp3:
          print "passed test: comma and notesglob."
          ipass += 1
     else:
          print "TEST FAIL!   comma and notesglob."
     # ================================================================== 
     if DEBUG:
          print
          print "  =>> OCULAR TEST, verify Display of tuple:"
          I.viewcomma( '*#tuple,test*', 'ytest', wild=False )
          print
     #    --------------------------------
     print "     Eyeball random string with queue POP=True:"
     print I.omaxcomma( "rand*,*ring", 'ytest', POP=True )
     print
     if I.lastkid( 'ytest' ) == lkid - 1 :
          print "passed test: POP deleted row as expected."
          ipass += 1
     else:
          print "TEST FAIL!   row not deleted per POP."
     lkid = I.lastkid( 'ytest' )
     print "     Current last kid PRIMARY KEY: ", lkid
     gotl = I.getkid( lkid, 'ytest' )
     if gotl == tmp2:
          print "passed test: getkid."
          ipass += 1
     else:
          print "TEST FAIL!   getkid."
     #    --------------------------------
     print "     Trying .select with default POP=False ..."
     got4 = I.select( '#tuple,test', 'ytest' )
     if got4[1] == 98:
          print "passed test: comma2list with wild=True."
          ipass += 1
     else:
          print "TEST FAIL!   comma2list with wild=True."
     # ================================================================== 
     print "     (Note: infile v0.50 has passed inspection.)"
     #  Test infile separately since it requires an external file.
     #  2009-09-14 v0.21:
     #     pzinfile working fine with text files.
     #     Round-trip on a binary file produced matching SHA256 signatures.
     print "     (Note: inweb  v0.50 has passed inspection.)"
     #  Test inweb separately since it requires an external website.
     #  2009-09-20 v0.22
     #     HTML from python.org appears fine with newlines preserved.
     print "----------------------------------------------------------------"
     print "DELETING rows older than 30 minutes from ytest."
     I.freshen( 0.0208, 'ytest' )
     #            ^= 30mins expressed in days.
     ipass += 1
     print "DROPPING table ytest2."
     I.droptable( 'ytest2' )
     ipass += 1
     #  print "     (Note: copysub and copycomma v0.50 have passed inspection.)"
     print "COPYING table ytest to ytest2."
     copy( '', 'ytest', 'ytest2' )
     #     however, they are not necessarily identical for kids may differ.
     ipass += 1
     print "     Assert copy and fifo methods:", 
     obj1 = True
     while obj1 != None:
          obj1 = I.fifo( 'ytest' )
          obj2 = I.fifo( 'ytest2' )
          assert obj1 == obj2, " :: tester: copy FAIL."
     print "PASSED, objects equivalent."
     ipass += 1
     assert I.lastkid('ytest') == 0, " :: tester: fifo FAIL."
     print "     (ytest and ytest2 should be empty due to iterated fifo.)"
     ipass += 1
     print "VACUUMing the entire database."
     I.clean()
     ipass += 1
     print "----------------------------------------------------------------"
     #  print "ipass =", ipass
     if ipass == 20:
          #      ^increment if you added a test ;-)
          print " *** tester    compiled: PASSED -- verify results above. ***"
     else:
          print " !!! tester     summary: FAILED! -- y_serial BROKEN."




def testfarm( dir=Farm.dir0, maxbarns=Farm.barns0, noobs=500 ):
     '''Test class Farm for bugs. Include path for directory.'''
     print "\n======================== testfarm =============================="
     ipass = 0
     if not DEBUG:
          print "[DEBUG switch is currently False.]"
     F = Farm( dir, maxbarns )
     print "     directory:", F.dir
     print "     maxbarns:", F.maxbarns
     ipass += 1
     door = 0
     F.farmin( 2009, 'testfarm farmin', 'ytest', door )
     D = Main( Farm.barn(F, door))
     assert 2009 == D.select(0, 'ytest'), "farmin: FAIL"
     ipass += 1
     testbarn = F.dir + 'testbarn.sqlite'
     F.reap( 'farmin', 'ytest', 'ytest', door, testbarn, wild=True )
     T = Main( testbarn )
     assert 2009 == T.select(0, 'ytest'), "reap: FAIL"
     ipass += 1
     print "passed: farmin and reap."
     D.droptable( 'ytest' )
     T.droptable( 'ytest' )
     print "----------------------------------------------------------------"
     print "TESTING plant:  %s fresh objects. [Stand-by ...]" % noobs
     for i in range( noobs ):
          F.plant( 'myobj', 'plant-%s' % i, 'ytest', testbarn )
     ipass += 1
     print "     lastkid in target database:", T.lastkid( 'ytest' )
     print "Next, reap_ALL_BARNS..."
     F.plant( 'reap_ALL_BARNS', '', 'ytest', testbarn )
     allkids = T.lastkid( 'ytest' )
     print "     lastkid in target database:", allkids
     assert noobs == allkids, "reap_ALL_BARNS:  FAIL."
     ipass += 1
     print "Cleaning up after plant objects."
     T.clean(     0, 'ytest' )
     F.cleanfarm( 0, 'ytest' )
     ipass += 1
     print "----------------------------------------------------------------"
     #  print "ipass =", ipass
     if ipass == 6:
          #      ^increment if you added a test ;-)
          print " *** testfarm  compiled: PASSED -- verify results above. ***"
     else:
          print " !!! testfarm   summary: FAILED! -- y_serial BROKEN."


if __name__ == "__main__":
     print "\n  ::  THIS IS A MODULE for import -- not for direct execution! \n"
     raw_input('Enter something to get out: ')


# ============================ Acknowledgements ====================================== 

#  - Special thanks to Bret TAYLOR for original inspiration and battle-tested 
#       case study of Friendfeed. 
#
#  - Alex MARTELLI for his assistance and his Cookbook for 
#       serving up delicious Python meals.


#  SQLite is a C library that provides a disk-based database that doesn't 
#  require a separate server process.  D. Richard HIPP deserves huge credit 
#  for putting SQLite in the public domain:  solving 80% of data persistance 
#  issues, using only 20% of the effort required by other SQL databases.
#  Thus SQLite is the most widely deployed SQL database engine in the world. 


#  The sqlite3 module was written by Gerhard HAERING and provides a SQL interface
#  compliant with the Python DB-API 2.0 specification described by `PEP 249
#  <http://www.python.org/dev/peps/pep-0249/>`_.  His effort -- updating that
#  module from version 2.3.2 in Python 2.5 to version 2.4.1 in Python 2.6 -- has
#  been crucial to y_serial. We are counting on him to resolve `bug 7723 
#  <http://bugs.python.org/issue7723>`_ for an easy transition to Python 3.


#  zlib was written by Jean-loup GAILLY (compression) and Mark ADLER
#  (decompression).  Jean-loup is also the primary author/maintainer of gzip(1);
#  Mark is also the author of gzip's and UnZip's main decompression routines and
#  was the original author of Zip.  Not surprisingly, the compression algorithm
#  used in zlib is essentially the same as that in gzip and Zip, namely, the
#  'deflate' method that originated in PKWARE's PKZIP 2.x.  


# =========================== Revised BSD License ==================================== 
#  
#  Copyright (c) 2009, y_Developers, http://yserial.sourceforge.net
#  All rights reserved.
#  
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#  
#      * Redistributions of source code must retain the above copyright notice, 
#        this list of conditions and the following disclaimer.
#  
#      * Redistributions in binary form must reproduce the above copyright notice, 
#        this list of conditions and the following disclaimer in the documentation 
#        and/or other materials provided with the distribution.
#  
#      * Neither the name of this organization nor the names of its 
#        contributors may be used to endorse or promote products derived from 
#        this software without specific prior written permission.
#  
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

