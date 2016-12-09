#!/usr/bin/python

from   integrator_utils.mysql import *
import os

class Node:
    def __init__ (self, id, node_type=None):
        self.id        = id
        self.parent    = None
        self.children  = []
        self.is_leaf   = False
        self.is_root   = False
        if node_type=="leaf":
            self.is_leaf   = True
        if node_type=="root":
            self.is_root   = True

    def build_tree (self, parent_of, children_of):
        if self.is_root:
            self.children = [Node(x) for x in children_of.keys() if not parent_of.has_key(x)]
        elif children_of.has_key(self.id):
            self.children = [Node(x) for x in children_of[self.id]]
        for child in self.children:
            child.parent = self
            child.build_tree (parent_of, children_of)

    def output(self, cursor=None, depth=0):
        name = ""
        if cursor:
            qry = "select displayName from Pathway where id = %s" % self.id
            rows = search_db(cursor, qry)
            if rows and rows[0]: name = rows[0][0]
        print "\t"*depth + str(self.id) + "  " + name
        for child in self.children:
            child.output(cursor, depth+1)
        return

#########################################
def reconstruct_pathway_tree(cursor):
    root = Node(-1, "root")
    qry =  "select PathwayHierarchy.pathwayId, PathwayHierarchy.childPathwayId "
    qry += "from PathwayHierarchy, Pathway  where "
    qry += "PathwayHierarchy.pathwayId=Pathway.id and Pathway.species = 'Homo sapiens'"
    rows = search_db(cursor,qry)
    parent_of = {}
    children_of = {}
    for row in rows:
        [parent,child] = row
        parent_of[child] = parent
        if not children_of.has_key(parent): children_of[parent] = []
        children_of[parent].append(child)

    root.build_tree (parent_of, children_of)

    return root


#########################################
def main():
    # note the skip-auto-rehash option in .ucsc_myql_conf
    # it is the equivalent to -A on the mysql command line
    # means: no autocompletion, which makes mysql get up mych faster
    db     = connect_to_mysql(user="cookiemonster", passwd=(os.environ['COOKIEMONSTER_PASSWORD']))
    if not db: exit(1)
    cursor = db.cursor()
    switch_to_db(cursor, 'reactome')

    pathway_tree = reconstruct_pathway_tree(cursor)
    pathway_tree.output(cursor)

    if False:
        qry  = "select id, stableId, displayName from Pathway where species='Homo sapiens'"
        rows = search_db(cursor,qry)
        for row in rows:
            print row
            [pathway_id, stableId, displayName] = row
            print displayName
            qry  = "select  ReactionLikeEvent_To_PhysicalEntity.physicalEntityId "
            qry += "from  Pathway_To_ReactionLikeEvent, ReactionLikeEvent_To_PhysicalEntity "
            qry += "where Pathway_To_ReactionLikeEvent.pathwayId=%s " % pathway_id
            qry += "and  Pathway_To_ReactionLikeEvent.reactionLikeEventId= ReactionLikeEvent_To_PhysicalEntity.reactionLikeEventId"
            rows2 = search_db(cursor,qry)
            physical_entity_ids = set([row2[0] for row2 in rows2])
            for physical_entity_id in physical_entity_ids:
                qry = "select displayName from PhysicalEntity where id=%s" % physical_entity_id
                entity_name = search_db(cursor,qry)[0][0]
                print "\t", physical_entity_id, entity_name

        exit(1)


    cursor.close()
    db.close()

    
    
    return True


#########################################
if __name__ == '__main__':
    main()


