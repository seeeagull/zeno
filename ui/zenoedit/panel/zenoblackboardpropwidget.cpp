#include "zenoblackboardpropwidget.h"
#include <zenoui/comctrl/zlabel.h>
#include <zenoui/style/zenostyle.h>
#include <QGridLayout>
#include "zenoapplication.h"
#include <zenoui/comctrl/zwidgetfactory.h>
#include <zenomodel/include/graphsmanagment.h>
#include "zenomainwindow.h"
#include <zenomodel/include/uihelper.h>

ZenoBlackboardPropWidget::ZenoBlackboardPropWidget(const QPersistentModelIndex &index, const QPersistentModelIndex &subIndex, QWidget *parent)
    : QWidget(parent), 
    m_idx(index), 
    m_subgIdx(subIndex) 
{
    QGridLayout *pGroupLayout = new QGridLayout(this);
    pGroupLayout->setContentsMargins(10, 15, 0, 15);
    pGroupLayout->setColumnStretch(1, 1);
    pGroupLayout->setColumnStretch(2, 3);
    pGroupLayout->setSpacing(10);
    PARAMS_INFO params = m_idx.data(ROLE_PARAMS_NO_DESC).value<PARAMS_INFO>();
    BLACKBOARD_INFO info = params["blackboard"].value.value<BLACKBOARD_INFO>();
    insertRow("background", PARAM_CONTROL::CONTROL_COLOR_NORMAL, info.background, 0, pGroupLayout);
    insertRow("title", PARAM_CONTROL::CONTROL_STRING, info.title, 1, pGroupLayout);
    insertRow("content", PARAM_CONTROL::CONTROL_MULTILINE_STRING, info.content, 2, pGroupLayout);
    IGraphsModel *pModel = zenoApp->graphsManagment()->currentModel();
    connect(pModel, SIGNAL(_dataChanged(const QModelIndex &, const QModelIndex &, int)), this,SLOT(onDataChanged(const QModelIndex &, const QModelIndex &, int)));
}

ZenoBlackboardPropWidget::~ZenoBlackboardPropWidget() 
{
}

void ZenoBlackboardPropWidget::onDataChanged(const QModelIndex &subGpIdx, const QModelIndex &idx, int role) {
    if (subGpIdx != m_subgIdx)
        return;
    if (role == ROLE_PARAMS_NO_DESC) {
        PARAMS_INFO params = idx.data(ROLE_PARAMS_NO_DESC).value<PARAMS_INFO>();
        BLACKBOARD_INFO info = params["blackboard"].value.value<BLACKBOARD_INFO>();
        m_pTextEdit->setText(info.content);
        m_pTitle->setText(info.title);
    }

}

void ZenoBlackboardPropWidget::insertRow(const QString &desc, const PARAM_CONTROL &ctrl, const QVariant &value, int row,QGridLayout *pGroupLayout) {    
    ZTextLabel *pLabel = new ZTextLabel(desc);
    pLabel->setFont(QFont("Segoe UI Semibold", 12));
    pLabel->setTextColor(QColor(255, 255, 255, 255 * 0.7));
    pLabel->setHoverCursor(Qt::ArrowCursor);

    ZIconLabel *pIcon = new ZIconLabel;
    pIcon->setIcons(ZenoStyle::dpiScaledSize(QSize(28, 28)), ":/icons/parameter_key-frame_idle.svg",
                    ":/icons/parameter_key-frame_hover.svg");
    pGroupLayout->addWidget(pIcon, row, 0, Qt::AlignCenter);

    pGroupLayout->addWidget(pLabel, row, 1, Qt::AlignLeft | Qt::AlignVCenter);

    CallbackCollection cbSet;
    cbSet.cbEditFinished = [=](QVariant newValue) {
        IGraphsModel *pModel = zenoApp->graphsManagment()->currentModel();
        if (!pModel)
            return;
        PARAMS_INFO params = m_idx.data(ROLE_PARAMS_NO_DESC).value<PARAMS_INFO>();
        BLACKBOARD_INFO info = params["blackboard"].value.value<BLACKBOARD_INFO>();
        if (desc == "title") {
            info.title = newValue.value<QString>();
        } else if (desc == "content") {
            info.content = newValue.value<QString>();
        } else if (desc == "background") {
            info.background = newValue.value<QColor>();
        }
        pModel->updateBlackboard(m_idx.data(ROLE_OBJID).toString(), info, m_subgIdx, true);
    };

    cbSet.cbSwitch = [=](bool bOn) {
        zenoApp->getMainWindow()->setInDlgEventLoop(bOn); //deal with ubuntu dialog slow problem when update viewport.
    };
    QWidget *pControl = zenoui::createWidget(value, ctrl, UiHelper::getControlDesc(ctrl), cbSet, QVariant()); 
    if (desc == "title" && dynamic_cast<ZLineEdit*>(pControl)) {
        m_pTitle = dynamic_cast<ZLineEdit *>(pControl);
    } else if (desc == "content" && dynamic_cast<ZTextEdit *>(pControl)) {
        m_pTextEdit = dynamic_cast<ZTextEdit *>(pControl);
    } 
    if (pControl)
        pGroupLayout->addWidget(pControl, row, 2, Qt::AlignVCenter);
}
