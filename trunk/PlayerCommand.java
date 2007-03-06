import java.io.*;

/**
 * A <code>PlayerCommand</code> is a <code>Command</code> that is send by a
 * client. This commandtype contains data about an action of a player,
 * there are several actions:
 * 1. move right;
 * 2. move left;
 * 3. move down;
 * 4. move up;
 * 5. nothing.
 * The nothing command is send at intervals, when a player has done nothing.
 */
public class PlayerCommand extends Command
{
	// Something
	public static final long serialVersionUID = 1;
	
	// The action id of the player's action
	private int action_id;
	
	/**
	 * Creates a new <code>PlayerCommand</code> object, with a
	 * given action and player id.
	 *
	 * @param sender_id The id of the player who send this command.
	 * @param action_id The id of the action.
	 */
	public PlayerCommand(int sender_id, int action_id)
	{
		super(Config.CMD_PLAYER, sender_id);
		this.action_id = action_id;
	}
	
	/**
	 * Returns the action id of this command.
	 *
	 * @return the action id of this <code>PlayerCommand</code>.
	 */
	public int getActionID()
	{
		return this.action_id;
	}
}