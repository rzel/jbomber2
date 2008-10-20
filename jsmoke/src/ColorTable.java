import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.table.*;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;

class ColorTable extends JTable {
	private static final long serialVersionUID = 1L; // prevent warning

	public ColorTable(DefaultTableModel m) {
		super(m);
	}

	public Class getColumnClass(int column) { //enable JTable to use different renderers, eg Checkbox for Boolean
		return getValueAt(0, column).getClass();
	}

	public void setPreferredColumnWidths(double[] percentages) {
		Dimension tableDim = this.getPreferredSize();
		double total = 0;
		for (int i = 0; i < getColumnModel().getColumnCount(); i++)
			total += percentages[i];
		for (int i = 0; i < getColumnModel().getColumnCount(); i++) {
			TableColumn column = getColumnModel().getColumn(i);
			column.setPreferredWidth((int)(tableDim.width * (percentages[i] / total)));
		}
	}

	public void setPreferredColumnWidths(int[] widths) {
		for (int i = 0; i < getColumnModel().getColumnCount(); i++) {
			TableColumn column = getColumnModel().getColumn(i);
			column.setPreferredWidth(widths[i]);
		}
	}
}


class ColorTableTableModelListener implements TableModelListener {
	public void tableChanged(TableModelEvent evt) {
		if (evt.getType() == TableModelEvent.UPDATE) {
			int column = evt.getColumn();
			int row = evt.getFirstRow();
		}
	}
}


class ColorTableRenderer extends DefaultTableCellRenderer { //for coloring cells
	private static final long serialVersionUID = 1L; // prevent warning

	public ColorTableRenderer() {
		super();
	}

	public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
		if (value != null) {
			setBackground((Color)value);
		} else {
			setText("Error!");
		}
		return this;
	}
}

class ColorTableRightClick extends MouseAdapter {
	private JMenuItem menuitemAddNewColor;
	private JMenuItem menuitemRemoveColor;
	private JMenuItem menuitemChooseColor;

	public ColorTableRightClick(ActionListener al) {
		menuitemAddNewColor = new JMenuItem("Add new color here");
		menuitemAddNewColor.addActionListener(al);
		menuitemAddNewColor.setActionCommand("COLORTABLE_ADD_NEW_COLOR");
		menuitemRemoveColor = new JMenuItem("Remove this color");
		menuitemRemoveColor.addActionListener(al);
		menuitemRemoveColor.setActionCommand("COLORTABLE_REMOVE_COLOR");
		menuitemChooseColor = new JMenuItem("Pick new color");
		menuitemChooseColor.addActionListener(al);
		menuitemChooseColor.setActionCommand("COLORTABLE_PICK_COLOR");
	}

	private void checkRightClick(MouseEvent e) {
		if (e.isPopupTrigger()) {
			JTable source = (JTable)e.getSource();
			int row = source.rowAtPoint(e.getPoint());
			int column = source.columnAtPoint(e.getPoint());
			source.changeSelection(row, column, false, false);
			JPopupMenu popup = new JPopupMenu();
			popup.add(menuitemAddNewColor);
			popup.add(menuitemRemoveColor);
			popup.add(menuitemChooseColor);
			menuitemRemoveColor.setEnabled(source.getRowCount() > 2); // NOTE: at least 2 colors for JSmoke!
			popup.show(e.getComponent(), e.getX(), e.getY());
		}
	}

	//Note: Popup menus are triggered differently on different systems.
	//Therefore, isPopupTrigger should be checked in both mousePressed and mouseReleased for proper cross-platform functionality.
	public void mouseReleased(MouseEvent e) { checkRightClick(e); }
	public void mousePressed(MouseEvent e) { checkRightClick(e); }
}