#ifndef __DESIGNER_COMMON_H__
#define __DESIGNER_COMMON_H__

enum NODE_ID
{
    HEADER,
    COMP_NODENAME,
    COMP_STATUS,
    COMP_CONTROL,
    COMP_HEADER_BACKBOARD,
    COMP_DISPLAY,

    BODY,
    COMP_LTSOCKET,
    COMP_LBSOCKET,
    COMP_RTSOCKET,
    COMP_RBSOCKET,
    COMP_BODYBACKBOARD,

    ELEMENT,
};

enum DESIGNER_NODE_TYPE
{
    NT_HIGHLAYER,
    NT_COMPONENT,
    NT_COMPONENT_AS_ELEMENT,
    NT_ELEMENT
};

enum NODE_CONTENT
{
    NC_NONE,
    NC_IMAGE,
    NC_BACKGROUND,
    NC_TEXT
};

enum NODEITEMROLE
{
    NODEID_ROLE = Qt::UserRole + 1,
    NODEPOS_ROLE,
    NODELOCK_ROLE,
    NODEVISIBLE_ROLE,

    NODEPATH_ROLE,
    NODEHOVERPATH_ROLE,
    NODESELECTEDPATH_ROLE,

    NODEFONT_ROLE,
    NODEFONTCOLOR_ROLE,

    NODETEXT_ROLE,
    NODETYPE_ROLE,
    NODECONTENT_ROLE,

    //BACKGROUND ROLE
    NODECOLOR_NORMAL_ROLE,
    NODECOLOR_HOVERD_ROLE,
    NODECOLOR_SELECTED_ROLE,

    NODE_LTRADIUS_ROLE,
    NODE_RTRADIUS_ROLE,
    NODE_LBRADIUS_ROLE,
    NODE_RBRADIUS_ROLE,
};

enum SnapWay {
    NO_SNAP,
    SNAP_GRID,
    SNAP_PIXEL
};

#define NODE_MODEL_NAME "nodemodeldata"
#define NODE_SELECTION_MODEL "nodeselectionmodel"

#define HEADER_ID "Header"
#define BODY_ID "Body"

#define PIXELS_IN_CELL 8

#define ZVALUE_GRID_BACKGROUND -100
#define ZVALUE_GRID_SMALL -11
#define ZVALUE_GRID_BIG -10

#define ZVALUE_CORE_ITEM -9

#define ZVALUE_LOCKED_BG -8
#define ZVALUE_LOCKED_CP -7
#define ZVALUE_LOCKED_ELEM -6

#define ZVALUE_BACKGROUND 0
#define ZVALUE_COMPONENT 5
#define ZVALUE_ELEMENT 10

#define ZVALUE_SELECTED 100

#endif