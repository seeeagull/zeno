#include "qdmgraphicsnode.h"
#include "qdmgraphicssocket.h"
#include "qdmgraphicsscene.h"
#include "qdmnodeparamedit.h"
#include <zeno/dop/Descriptor.h>
#include <zeno/ztd/algorithm.h>
#include <zeno/zmt/log.h>
#include <QLineEdit>
#include <algorithm>

ZENO_NAMESPACE_BEGIN

QDMGraphicsNode::QDMGraphicsNode()
{
    setFlag(QGraphicsItem::ItemIsMovable);
    setFlag(QGraphicsItem::ItemIsSelectable);

    label = new QGraphicsTextItem(this);
    label->setDefaultTextColor(QColor(0xcccccc));
    label->setPos(0, -SOCKSTRIDE);
}

QDMGraphicsNode::~QDMGraphicsNode() = default;

float QDMGraphicsNode::getHeight() const
{
    size_t count = std::max(socketIns.size(), socketOuts.size());
    return SOCKMARGINTOP + std::max(SOCKSTRIDE * count, MINHEIGHT) + SOCKMARGINBOT;
}

QRectF QDMGraphicsNode::boundingRect() const
{
    return {-QDMGraphicsSocket::SIZE, 1e-6f, WIDTH + QDMGraphicsSocket::SIZE * 2, getHeight() - 2e-6f};
}

void QDMGraphicsNode::paint(QPainter *painter, QStyleOptionGraphicsItem const *styleOptions, QWidget *widget)
{
    if (isSelected()) {
        QPen pen;
        pen.setColor(QColor(0xff8800));
        pen.setWidthF(BORDER);
        painter->setPen(pen);
    } else {
        painter->setPen(Qt::NoPen);
    }
    painter->setBrush(QColor(0x555555));

    QPainterPath path;
    QRectF rect(0, 0, WIDTH, getHeight());
    path.addRoundedRect(rect, ROUND, ROUND);
    painter->drawPath(path.simplified());
}

QDMGraphicsSocketIn *QDMGraphicsNode::addSocketIn()
{
    auto socketIn = new QDMGraphicsSocketIn;
    socketIn->setParentItem(this);

    size_t index = socketIns.size();
    socketIn->setPos(-socketIn->SIZE / 2, SOCKMARGINTOP + SOCKSTRIDE * index);

    socketIns.push_back(socketIn);
    return socketIn;
}

QDMGraphicsSocketOut *QDMGraphicsNode::addSocketOut()
{
    auto socketOut = new QDMGraphicsSocketOut;
    socketOut->setParentItem(this);

    size_t index = socketOuts.size();
    socketOut->setPos(WIDTH + socketOut->SIZE / 2, SOCKMARGINTOP + SOCKSTRIDE * index);

    socketOuts.push_back(socketOut);
    return socketOut;
}

void QDMGraphicsNode::initByType(const std::string &type)
{
    auto const &desc = dop::descriptor_table().at(type);
    initByDescriptor(&desc);
}

void QDMGraphicsNode::initByDescriptor(const dop::Descriptor *desc)
{
    this->desc = desc;
    for (auto const &sockinfo: desc->inputs) {
        auto socket = addSocketIn();
        socket->setName(sockinfo.name);
        socket->setType(sockinfo.type);
        socket->setDefl(sockinfo.defl);
    }
    for (auto const &sockinfo: desc->outputs) {
        auto socket = addSocketOut();
        socket->setName(sockinfo.name);
        socket->setType(sockinfo.type);
        socket->setDefl(sockinfo.defl);
    }

    auto name = getScene()->allocateNodeName(desc->name);
    setName(QString::fromStdString(name));
}

void QDMGraphicsNode::setName(QString name)
{
    label->setPlainText(name);
    this->name = name.toStdString();
}

QDMGraphicsSocketIn *QDMGraphicsNode::socketInAt(size_t index)
{
    return socketIns.at(index);
}

QDMGraphicsSocketOut *QDMGraphicsNode::socketOutAt(size_t index)
{
    return socketOuts.at(index);
}

QDMGraphicsScene *QDMGraphicsNode::getScene() const
{
    return static_cast<QDMGraphicsScene *>(scene());
}

QDMGraphicsScene *QDMGraphicsNode::getSubnetScene() const
{
    return nullptr;
}

void QDMGraphicsNode::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
    if (event->button() == Qt::RightButton) {
        getScene()->removeNode(this);
        return;
    }

    QGraphicsItem::mousePressEvent(event);
}

void QDMGraphicsNode::unlinkAll()
{
    for (auto const &p: socketIns) {
        p->unlinkAll();
    }
    for (auto const &p: socketOuts) {
        p->unlinkAll();
    }
}

size_t QDMGraphicsNode::socketInIndex(QDMGraphicsSocketIn *socket)
{
    return ztd::find_index(socketIns, socket);
}

size_t QDMGraphicsNode::socketOutIndex(QDMGraphicsSocketOut *socket)
{
    return ztd::find_index(socketOuts, socket);
}

/*void QDMGraphicsNode::socketUnlinked(QDMGraphicsSocketIn *socket)
{
    auto &sockIn = dopNode->inputs.at(socketInIndex(socket));
    sockIn.node = nullptr;
    sockIn.sockid = 0;
}

void QDMGraphicsNode::socketLinked(QDMGraphicsSocketIn *socket, QDMGraphicsSocketOut *srcSocket)
{
    auto srcNode = static_cast<QDMGraphicsNode *>(srcSocket->parentItem());
    auto &sockIn = dopNode->inputs.at(socketInIndex(socket));
    sockIn.node = srcNode->dopNode.get();
    sockIn.sockid = srcNode->socketOutIndex(srcSocket);
}

void QDMGraphicsNode::socketValueChanged(QDMGraphicsSocketIn *socket)
{
    auto &sockIn = dopNode->inputs.at(socketInIndex(socket));
    sockIn.value = socket->value;
}*/

void QDMGraphicsNode::invalidate()
{
    ZENO_INFO("invalidate node");
}

void QDMGraphicsNode::setupParamEdit(QDMNodeParamEdit *paredit)
{
    for (size_t i = 0; i < socketIns.size(); i++) {
        auto sockIn = socketIns[i];
        if (auto edit = paredit->makeEditForType(this, i, sockIn->getType()))
            paredit->addRow(QString::fromStdString(sockIn->getName()), edit);
    }
}

std::vector<std::string> QDMGraphicsNode::getInputNames() const {
    std::vector<std::string> ret;
    for (auto const &in: socketIns)
        ret.push_back(in->getName());
    return ret;
}

std::vector<std::string> QDMGraphicsNode::getOutputNames() const {
    std::vector<std::string> ret;
    for (auto const &out: socketOuts)
        ret.push_back(out->getName());
    return ret;
}

QDMGraphicsNode *QDMGraphicsNode::underlyingNode() {
    return this;
}

ZENO_NAMESPACE_END
