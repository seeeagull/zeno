#include <QTextEdit>
#include "zfxhighlighter.h"
#include "zfxkeywords.h"

ZfxHighlighter::ZfxHighlighter(QTextEdit* textEdit)
	: QSyntaxHighlighter(textEdit)
	, m_pTextEdit(textEdit)
{
	initRules();
	if (!m_pTextEdit) return;
	m_pTextEdit->installEventFilter(this);
	connect(m_pTextEdit, &QTextEdit::cursorPositionChanged, this, &ZfxHighlighter::highlightCurrentLine);
	connect(m_pTextEdit, &QTextEdit::selectionChanged, this, &ZfxHighlighter::onSelectionChanged);
}

void ZfxHighlighter::highlightBlock(const QString& text)
{
	if (text.size() == 0) return;
	for (const HighlightingRule& rule : m_highlightingRules) {
		QRegularExpressionMatchIterator matchIterator = rule.pattern.globalMatch(text);
		while (matchIterator.hasNext()) {
			QRegularExpressionMatch match = matchIterator.next();
			setFormat(match.capturedStart(), match.capturedLength(), rule.format);
		}
	}
}

void ZfxHighlighter::initRules()
{
	HighlightingRule rule;

	// number
	rule.pattern = QRegularExpression("\\b\\d+(\\.\\d+)?\\b");
	rule.format = m_highlightTheme.format(ZfxTextStyle::C_NUMBER);
	m_highlightingRules.append(rule);

	// string, seems no need

	// local var
	rule.pattern = QRegularExpression("\\@[a-zA-Z]+\\b");
	rule.format = m_highlightTheme.format(ZfxTextStyle::C_LOCALVAR);
	m_highlightingRules.append(rule);

	// global var
	rule.pattern = QRegularExpression("\\$[a-zA-Z]+\\b");
	rule.format = m_highlightTheme.format(ZfxTextStyle::C_GLOBALVAR);
	m_highlightingRules.append(rule);

	// keyword

	// func, only highlight supported
	QStringList functionPatterns;
	for (const auto& func : zfxFunction) {
		functionPatterns << QString("\\b%1\\b").arg(func);
	}
	auto functionFormat = m_highlightTheme.format(ZfxTextStyle::C_FUNCTION);
	for (const auto& pattern : functionPatterns) {
		rule.pattern = QRegularExpression(pattern);
		rule.format = functionFormat;
		m_highlightingRules.append(rule);
	}

	// comment
	rule.pattern = QRegularExpression("#[^\n]*");
	rule.format = m_highlightTheme.format(ZfxTextStyle::C_COMMENT);
	m_highlightingRules.append(rule);
}

void ZfxHighlighter::highlightCurrentLine()
{
	if (!m_pTextEdit || m_pTextEdit->isReadOnly() || !m_pTextEdit->hasFocus()) return;

	QList<QTextEdit::ExtraSelection> extraSelections;
	QTextEdit::ExtraSelection selection;
	selection.format = m_highlightTheme.format(ZfxTextStyle::C_CurrentLine);
	selection.format.setProperty(QTextFormat::FullWidthSelection, true);
	selection.cursor = m_pTextEdit->textCursor();
	selection.cursor.clearSelection();
	extraSelections.append(selection);
	m_pTextEdit->setExtraSelections(extraSelections);
}

void ZfxHighlighter::onSelectionChanged()
{
}

bool ZfxHighlighter::eventFilter(QObject* object, QEvent* event)
{
	if (object == m_pTextEdit) {
		if (event->type() == QEvent::FocusOut) {
			// clear selection
			QTextCursor cursor = m_pTextEdit->textCursor();
			cursor.clearSelection();
			m_pTextEdit->setTextCursor(cursor);
			m_pTextEdit->setExtraSelections(QList<QTextEdit::ExtraSelection>());
		}
		else if (event->type() == QEvent::FocusIn) {
			highlightCurrentLine();
		}
	}
	return QSyntaxHighlighter::eventFilter(object, event);
}
