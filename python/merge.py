import sqlite3, sys
#
## usage python merge.py src dst
#
con = sqlite3.connect(sys.argv[2])
cur = con.cursor()
sql = "ATTACH '"+sys.argv[1]+"' as src"
cur.execute(sql)
cur.close()
cur = con.cursor()
sql = "SELECT * FROM main.sqlite_master WHERE type='table'"
cur.execute(sql)
main_tables = cur.fetchall()
cur.close()
cur = con.cursor()
sql = "SELECT * FROM src.sqlite_master WHERE type='table'"
cur.execute(sql)
src_tables = cur.fetchall()
cur.close()
for var in src_tables:
  varname = var[1]
  if varname not in [x[1] for x in src_tables]:
    cur = con.cursor()
    cur.execute(var[-1])
    cur.close()
  cur = con.cursor()
  sql = "INSERT OR REPLACE into "+varname+" SELECT * FROM src."+varname
  cur.execute(sql)
  cur.close()
con.commit()
con.close()
exit(0)
