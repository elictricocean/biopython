

init_sql = """
DROP TABLE IF EXISTS term;

CREATE TABLE term (
  id int(11) NOT NULL,
  name varchar(255) NOT NULL default '',
  term_type varchar(55) NOT NULL,
  acc varchar(255) NOT NULL,
  is_obsolete int(11) NOT NULL default '0',
  is_root int(11) NOT NULL default '0',
  is_relation int(11) NOT NULL default '0',
  PRIMARY KEY (id) 
);

CREATE UNIQUE INDEX acc ON term (acc);
CREATE INDEX t1 ON term (name);
CREATE INDEX t2 ON term (term_type);
CREATE INDEX t3 ON term (acc);
CREATE INDEX t4 ON term (id, acc);
CREATE INDEX t5 ON term (id, name);
CREATE INDEX t6 ON term (id, term_type);
CREATE INDEX t7 ON term (id, acc, name, term_type);

DROP TABLE IF EXISTS term2term;
CREATE TABLE term2term (
  id int(11),
  relationship_type_id int(11) NOT NULL,
  term1_id int(11) NOT NULL,
  term2_id int(11) NOT NULL,
  complete int(11) NOT NULL default '0',
  PRIMARY KEY  (id)
);
CREATE UNIQUE INDEX term1_id ON term2term (term1_id, term2_id, relationship_type_id);
CREATE INDEX tt1 ON term2term (term1_id);
CREATE INDEX tt2 ON term2term (term2_id);
CREATE INDEX tt3 ON term2term (term1_id, term2_id);
CREATE INDEX tt4 ON term2term (relationship_type_id);

DROP TABLE IF EXISTS term_synonym;
CREATE TABLE term_synonym (
  term_id int(11) NOT NULL,
  term_synonym varchar(330),
  acc_synonym varchar(255),
  synonym_type_id int(11) NOT NULL,
  synonym_category_id int(11)
);
CREATE UNIQUE INDEX term_id ON term_synonym (term_id, term_synonym);
CREATE INDEX synonym_type_id ON term_synonym (synonym_type_id);
CREATE INDEX synonym_category_id ON term_synonym (synonym_category_id);
CREATE INDEX ts1 ON term_synonym (term_id);
CREATE INDEX ts2 ON term_synonym (term_synonym);
CREATE INDEX ts3 ON term_synonym (term_id, term_synonym);
"""


class OntologySQLWriter(object):
    
    def __init__(self, conn):
        self.conn = conn
        self.cursor = self.conn.cursor()
        for sql in init_sql.split(";"):
            self.cursor.execute(sql)
        self.id_map = { 'is_a' : 0 }
        insert_term = "INSERT OR REPLACE INTO term(id, name, acc, term_type) VALUES(0,'is_a','relationship','is_a')"
        self.cursor.execute(insert_term)
    
    def write_record(self, record):
        if record.tags['id'][0] not in self.id_map:
            self.id_map[record.tags['id'][0]] = len(self.id_map)
        insert_term = "INSERT OR REPLACE INTO term(id, name, acc, term_type) VALUES(?,?,?,?)"
        term_type = ""
        if record.name == "Term":
            term_type = "gene_ontology"
        if record.name == "Typedef":
            term_type = "relationship"
        
        self.cursor.execute(insert_term, 
            [ 
                self.id_map[record.tags['id'][0]], 
                str(record.tags['name'][0]), 
                str(record.tags['id'][0]),
                term_type
            ] 
        )
        
        for tag in record.tags:
            if tag == "is_a":
                for a in record.tags["is_a"]:
                    if a not in self.id_map:
                        self.id_map[a] = len(self.id_map)                        
                    insert_term2term = "INSERT INTO term2term(relationship_type_id, term1_id, term2_id) VALUES(?,?,?)"
                    self.cursor.execute(insert_term2term,
                        [
                            self.id_map['is_a'], 
                            self.id_map[record.tags['id'][0]], 
                            self.id_map[a]
                        ]
                    )
                
        self.conn.commit()
    
